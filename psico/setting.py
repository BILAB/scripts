'''
(c) 2011 Thomas Holder, MPI for Developmental Biology

License: BSD-2-Clause
'''

from pymol import cmd

def save_settings(filename='~/.pymolrc-settings.py', quiet=1):
    '''
DESCRIPTION
 
    Dumps all settings with non-default values to ~/.pymolrc-settings.py
 
    Feature Request: Save settings for later use - ID: 1009951
    https://sourceforge.net/tracker/?func=detail&aid=1009951&group_id=4546&atid=354546
    '''
    from pymol.setting import get_name_list
    quiet = int(quiet)
    if not filename.endswith('.py'):
        print('Warning: filename should end with ".py"')
    # temporatily load default settings and remember them
    cmd.reinitialize('store_defaults')
    cmd.reinitialize('original_settings')
    original = [(name, cmd.get(name)) for name in get_name_list()]
    cmd.reinitialize('settings')
    # dump to file
    filename = cmd.exp_path(filename)
    f = open(filename, 'w')
    f.write('# AUTOGENERATED FILE\n')
    f.write('from pymol import cmd, invocation\n')
    f.write('if invocation.options.show_splash: ') # no newline
    f.write('    print("Loading settings from " + ' + repr(filename) + ')\n')
    count = 0
    for name, o_value in original:
        value = cmd.get(name)
        if value != o_value:
            f.write('cmd.set("%s", %s)\n' % (name, repr(value)))
            if not quiet:
                print('set %s, %s # default: %s' % (name, value, o_value))
            count += 1
    f.close()
    if not quiet:
        print('Dumped %d settings to %s' % (count, filename))

def paper_settings(fancy=0, quiet=0):
    '''
DESCRIPTION

    Set rendering quality high and some stuff good for printing:
    
     * Side chain helper (cartoon_side_chain_helper = 1)
     * No shadows (ray_shadows = 0)

ARGUMENTS

    fancy = 0 or 1: set cartoon_fancy_helices and cartoon_highlight_color
    {default: 0}

NOTES

    You may also try "set ray_trace_mode, 1"
    '''
    fancy, quiet = int(fancy), int(quiet)
    if fancy == 1:
        cmd.set('cartoon_fancy_helices', 1, quiet=quiet)
        cmd.set('cartoon_highlight_color', 'grey50', quiet=quiet)
    cmd.set('cartoon_side_chain_helper', 1, quiet=quiet)
    cmd.set('ray_shadows', 0, quiet=quiet)
    cmd.set('opaque_background', 0, quiet=quiet)
    cmd.bg_color('white')

class set_temporary(object):
    '''
DESCRIPTION

    API only. Supports the following pattern:

    >>> with set_temporary(pdb_retain_ids=1):
    ...    cmd.save('out.pdb')
    '''
    def __init__(self, *args, **kwargs):
        self.sele = kwargs.pop('selection', '')
        self.args = args + tuple(kwargs.items())
    def __enter__(self):
        self.saved = []
        for k, v in self.args:
            v_saved = cmd.get(k, self.sele)
            if v != v_saved:
                self.saved.append((k, v_saved))
                cmd.set(k, v, self.sele)
        return self
    def __exit__(self, type, value, traceback):
        for k, v in self.saved:
            cmd.set(k, v, self.sele)

cmd.extend('save_settings', save_settings)
cmd.extend('paper_settings', paper_settings)

# vi:expandtab:smarttab