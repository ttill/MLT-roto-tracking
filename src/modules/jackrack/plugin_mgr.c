/*
 *   jack-ladspa-host
 *    
 *   Copyright (C) Robert Ham 2002, 2003 (node@users.sourceforge.net)
 *    
 *   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program; if not, write to the Free Software
 *   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <dlfcn.h>
#include <math.h>
#include <strings.h>
#include <ctype.h>
#include <ladspa.h>

#include "plugin_mgr.h"
#include "plugin_desc.h"
#include "ui.h"


static gboolean
plugin_is_valid (const LADSPA_Descriptor * descriptor)
{
  unsigned long i;
  unsigned long icount = 0;
  unsigned long ocount = 0;
  
  for (i = 0; i < descriptor->PortCount; i++)
    {
      if (!LADSPA_IS_PORT_AUDIO (descriptor->PortDescriptors[i]))
        continue;
      
      if (LADSPA_IS_PORT_INPUT (descriptor->PortDescriptors[i]))
        icount++;
      else
        ocount++;
    }
  
  if (icount == 0 || ocount == 0)
    return FALSE;
  
  return TRUE;
}

static void
plugin_mgr_get_object_file_plugins (ui_t * ui, plugin_mgr_t * plugin_mgr, const char * filename)
{
  const char * dlerr;
  void * dl_handle;
  LADSPA_Descriptor_Function get_descriptor;
  const LADSPA_Descriptor * descriptor;
  unsigned long plugin_index;
  plugin_desc_t * desc, * other_desc = NULL;
  GSList * list;
  gboolean exists;
  int err;
  
  /* open the object file */
  dl_handle = dlopen (filename, RTLD_NOW|RTLD_GLOBAL);
  if (!dl_handle)
    {
      fprintf (stderr, "%s: error opening shared object file '%s': %s\n",
               __FUNCTION__, filename, dlerror());
      return;
    }
  
  
  /* get the get_descriptor function */
  dlerror (); /* clear the error report */
  
  get_descriptor = (LADSPA_Descriptor_Function)
    dlsym (dl_handle, "ladspa_descriptor");
  
  dlerr = dlerror();
  if (dlerr) {
    fprintf (stderr, "%s: error finding ladspa_descriptor symbol in object file '%s': %s\n",
             __FUNCTION__, filename, dlerr);
    dlclose (dl_handle);
    return;
  }
  
  plugin_index = 0;
  while ( (descriptor = get_descriptor (plugin_index)) )
    {
      if (!plugin_is_valid (descriptor))
        {
          plugin_index++;
          continue;
        }

      
      /* check it doesn't already exist */
      exists = FALSE;
      for (list = plugin_mgr->all_plugins; list; list = g_slist_next (list))
        {
          other_desc = (plugin_desc_t *) list->data;
          
          if (other_desc->id == descriptor->UniqueID)
            {
              exists = TRUE;
              break;
            }
        }
      
      if (exists)
        {
          printf ("Plugin %ld exists in both '%s' and '%s'; using version in '%s'\n",
                  descriptor->UniqueID, other_desc->object_file, filename, other_desc->object_file);
          plugin_index++;
          continue;
        }

      
      desc = plugin_desc_new_with_descriptor (filename, plugin_index, descriptor);
      plugin_mgr->all_plugins = g_slist_append (plugin_mgr->all_plugins, desc);
      plugin_index++;
      plugin_mgr->plugin_count++;
      
      /* print in the splash screen */
      /* printf ("Loaded plugin '%s'\n", desc->name); */
    }
  
  err = dlclose (dl_handle);
  if (err)
    {
      fprintf (stderr, "%s: error closing object file '%s': %s\n",
               __FUNCTION__, filename, dlerror ());
    }
}

static void
plugin_mgr_get_dir_plugins (ui_t * ui, plugin_mgr_t * plugin_mgr, const char * dir)
{
  DIR * dir_stream;
  struct dirent * dir_entry;
  char * file_name;
  int err;
  size_t dirlen;
  
  dir_stream = opendir (dir);
  if (!dir_stream)
    {
/*      fprintf (stderr, "%s: error opening directory '%s': %s\n",
               __FUNCTION__, dir, strerror (errno)); */
      return;
    }
  
  dirlen = strlen (dir);
  
  while ( (dir_entry = readdir (dir_stream)) )
    {
      if (strcmp (dir_entry->d_name, ".") == 0 ||
          strcmp (dir_entry->d_name, "..") == 0)
        continue;
  
      file_name = g_malloc (dirlen + 1 + strlen (dir_entry->d_name) + 1);
    
      strcpy (file_name, dir);
      if (file_name[dirlen - 1] == '/')
        strcpy (file_name + dirlen, dir_entry->d_name);
      else
        {
          file_name[dirlen] = '/';
          strcpy (file_name + dirlen + 1, dir_entry->d_name);
        }
    
      plugin_mgr_get_object_file_plugins (ui, plugin_mgr, file_name);
      
      g_free (file_name);
    }

  err = closedir (dir_stream);
  if (err)
    fprintf (stderr, "%s: error closing directory '%s': %s\n",
             __FUNCTION__, dir, strerror (errno));
}

static void
plugin_mgr_get_path_plugins (ui_t * ui, plugin_mgr_t * plugin_mgr)
{
  char * ladspa_path, * dir;
  
  ladspa_path = g_strdup (getenv ("LADSPA_PATH"));
  if (!ladspa_path)
    ladspa_path = g_strdup ("/usr/local/lib/ladspa:/usr/lib/ladspa");
  
  dir = strtok (ladspa_path, ":");
  do
    plugin_mgr_get_dir_plugins (ui, plugin_mgr, dir);
  while ((dir = strtok (NULL, ":")));

  g_free (ladspa_path);
}

static gint
plugin_mgr_sort (gconstpointer a, gconstpointer b)
{
  const plugin_desc_t * da;
  const plugin_desc_t * db;
  da = (const plugin_desc_t *) a;
  db = (const plugin_desc_t *) b;
  
  return strcasecmp (da->name, db->name);
}

plugin_mgr_t *
plugin_mgr_new (ui_t * ui)
{
  plugin_mgr_t * pm;
  
  pm = g_malloc (sizeof (plugin_mgr_t));
  pm->all_plugins = NULL;  
  pm->plugins = NULL;
  pm->plugin_count = 0;
  
  plugin_mgr_get_path_plugins (ui, pm);
  
  if (!pm->all_plugins)
    {
      fprintf (stderr, "No LADSPA plugins were found!\n\nCheck your LADSPA_PATH environment variable.\n");
      abort ();
    }
  
  pm->all_plugins = g_slist_sort (pm->all_plugins, plugin_mgr_sort);
  
  return pm;
}

void
plugin_mgr_destroy (plugin_mgr_t * plugin_mgr)
{
  GSList * list;
  
  for (list = plugin_mgr->all_plugins; list; list = g_slist_next (list))
    plugin_desc_destroy ((plugin_desc_t *) list->data);
  
  g_slist_free (plugin_mgr->plugins);
  g_slist_free (plugin_mgr->all_plugins);
  free (plugin_mgr);
}


void
plugin_mgr_set_plugins (plugin_mgr_t * plugin_mgr, unsigned long rack_channels)
{
  GSList * list;
  plugin_desc_t * desc;

  /* clear the current plugins */
  g_slist_free (plugin_mgr->plugins);
  plugin_mgr->plugins = NULL;
  
  for (list = plugin_mgr->all_plugins; list; list = g_slist_next (list))
    {
      desc = (plugin_desc_t *) list->data;
      
      if (plugin_desc_get_copies (desc, rack_channels) != 0)
        plugin_mgr->plugins = g_slist_append (plugin_mgr->plugins, desc);
    }
}

static plugin_desc_t *
plugin_mgr_find_desc (plugin_mgr_t * plugin_mgr, GSList * plugins, unsigned long id)
{
  GSList * list;
  plugin_desc_t * desc;
  
  for (list = plugins; list; list = g_slist_next (list))
    {
      desc = (plugin_desc_t *) list->data;
      
      if (desc->id == id)
        return desc;
    }
  
  return NULL;
}

plugin_desc_t *
plugin_mgr_get_desc (plugin_mgr_t * plugin_mgr, unsigned long id)
{
  return plugin_mgr_find_desc (plugin_mgr, plugin_mgr->plugins, id);
}

plugin_desc_t *
plugin_mgr_get_any_desc (plugin_mgr_t * plugin_mgr, unsigned long id)
{
  return plugin_mgr_find_desc (plugin_mgr, plugin_mgr->all_plugins, id);
}


/* EOF */