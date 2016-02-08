#!/usr/bin/env python
'''
$Id: showMIweb.py 1531 2014-02-05 15:02:36Z apandini $

Copyright (C) 2009-13  Alessandro Pandini

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import sys
import os
import fnmatch
import time
import Tkinter
from Tkinter import *
import Pmw
from pymol import cmd
from pymol import stored

def __init__(self):
    self.menuBar.addmenuitem(   'Plugin', 'command',
                                'showMIweb',
                                label = 'showMIweb',
                                command = lambda s=self : DistMatTool(s))

class PmwFileDialog(Pmw.Dialog):
    """File Dialog using Pmw"""
    def __init__(self, parent = None, **kw):
	# Define the megawidget options.
	optiondefs = (
	    ('filter',    '*',              self.newfilter),
	    ('directory', os.getcwd(),      self.newdir),
	    ('filename',  '',               self.newfilename),
	    ('historylen',10,               None),
	    ('command',   None,             None),
            ('info',      None,             None),
	    )
	self.defineoptions(kw, optiondefs)
        # Initialise base class (after defining options).
	Pmw.Dialog.__init__(self, parent)

	self.withdraw()

        # Create the components.
	interior = self.interior()

        if self['info'] is not None:
            rowoffset=1
            dn = self.infotxt()
            dn.grid(row=0,column=0,columnspan=2,padx=3,pady=3)
        else:
            rowoffset=0

	dn = self.mkdn()
	dn.grid(row=0+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	del dn

	# Create the directory list component.
	dnb = self.mkdnb()
	dnb.grid(row=1+rowoffset,column=0,sticky='news',padx=3,pady=3)
	del dnb

	# Create the filename list component.
	fnb = self.mkfnb()
	fnb.grid(row=1+rowoffset,column=1,sticky='news',padx=3,pady=3)
	del fnb

	# Create the filter entry
	ft = self.mkft()
	ft.grid(row=2+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	del ft

	# Create the filename entry
	fn = self.mkfn()
	fn.grid(row=3+rowoffset,column=0,columnspan=2,padx=3,pady=3)
	fn.bind('<Return>',self.okbutton)
	del fn

	# Buttonbox already exists
	bb=self.component('buttonbox')
	bb.add('OK',command=self.okbutton)
	bb.add('Cancel',command=self.cancelbutton)
	del bb

	Pmw.alignlabels([self.component('filename'),
			 self.component('filter'),
			 self.component('dirname')])

    def infotxt(self):
        """ Make information block component at the top """
        return self.createcomponent(
                'infobox',
                (), None,
                Tkinter.Label, (self.interior(),),
                width=51,
                relief='groove',
                foreground='darkblue',
                justify='left',
                text=self['info']
            )

    def mkdn(self):
        """Make directory name component"""
        return self.createcomponent(
	    'dirname',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['directory'],
	    entryfield_entry_width=40,
            entryfield_validate=self.dirvalidate,
	    selectioncommand=self.setdir,
	    labelpos='w',
	    label_text='Directory:')

    def mkdnb(self):
        """Make directory name box"""
        return self.createcomponent(
	    'dirnamebox',
	    (), None,
	    Pmw.ScrolledListBox, (self.interior(),),
	    label_text='directories',
	    labelpos='n',
	    hscrollmode='none',
	    dblclickcommand=self.selectdir)

    def mkft(self):
        """Make filter"""
        return self.createcomponent(
	    'filter',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['filter'],
	    entryfield_entry_width=40,
	    selectioncommand=self.setfilter,
	    labelpos='w',
	    label_text='Filter:')

    def mkfnb(self):
        """Make filename list box"""
        return self.createcomponent(
	    'filenamebox',
	    (), None,
	    Pmw.ScrolledListBox, (self.interior(),),
	    label_text='files',
	    labelpos='n',
	    hscrollmode='none',
	    selectioncommand=self.singleselectfile,
	    dblclickcommand=self.selectfile)

    def mkfn(self):
        """Make file name entry"""
        return self.createcomponent(
	    'filename',
	    (), None,
	    Pmw.ComboBox, (self.interior(),),
	    entryfield_value=self['filename'],
	    entryfield_entry_width=40,
            entryfield_validate=self.filevalidate,
	    selectioncommand=self.setfilename,
	    labelpos='w',
	    label_text='Filename:')
    
    def dirvalidate(self,string):
        if os.path.isdir(string):
            return Pmw.OK
        else:
            return Pmw.PARTIAL
        
    def filevalidate(self,string):
        if string=='':
            return Pmw.PARTIAL
        elif os.path.isfile(string):
            return Pmw.OK
        elif os.path.exists(string):
            return Pmw.PARTIAL
        else:
            return Pmw.OK
        
    def okbutton(self):
	"""OK action: user thinks he has input valid data and wants to
           proceed. This is also called by <Return> in the filename entry"""
	fn=self.component('filename').get()
	self.setfilename(fn)
	if self.validate(fn):
	    self.canceled=0
	    self.deactivate()

    def cancelbutton(self):
	"""Cancel the operation"""
	self.canceled=1
	self.deactivate()

    def tidy(self,w,v):
	"""Insert text v into the entry and at the top of the list of 
           the combobox w, remove duplicates"""
	if not v:
	    return
	entry=w.component('entry')
	entry.delete(0,'end')
	entry.insert(0,v)
	list=w.component('scrolledlist')
	list.insert(0,v)
	index=1
	while index<list.index('end'):
	    k=list.get(index)
	    if k==v or index>self['historylen']:
		list.delete(index)
	    else:
		index=index+1
        w.checkentry()

    def setfilename(self,value):
	if not value:
	    return
	value=os.path.join(self['directory'],value)
	dir,fil=os.path.split(value)
	self.configure(directory=dir,filename=value)
        
	c=self['command']
	if callable(c):
	    c()

    def newfilename(self):
	"""Make sure a newly set filename makes it into the combobox list"""
	self.tidy(self.component('filename'),self['filename'])
	
    def setfilter(self,value):
	self.configure(filter=value)

    def newfilter(self):
	"""Make sure a newly set filter makes it into the combobox list"""
	self.tidy(self.component('filter'),self['filter'])
	self.fillit()

    def setdir(self,value):
	self.configure(directory=value)

    def newdir(self):
	"""Make sure a newly set dirname makes it into the combobox list"""
	self.tidy(self.component('dirname'),self['directory'])
	self.fillit()

    def singleselectfile(self):
	"""Single click in file listbox. Move file to "filename" combobox"""
	cs=self.component('filenamebox').curselection()
	if cs!=():
	    value=self.component('filenamebox').get(cs)
            self.setfilename(value)

    def selectfile(self):
	"""Take the selected file from the filename, normalize it, and OK"""
        self.singleselectfile()
	value=self.component('filename').get()
        self.setfilename(value)
        if value:
	    self.okbutton()

    def selectdir(self):
	"""Take selected directory from the dirnamebox into the dirname"""
	cs=self.component('dirnamebox').curselection()
	if cs!=():
	    value=self.component('dirnamebox').get(cs)
	    dir=self['directory']
	    if not dir:
		dir=os.getcwd()
	    if value:
		if value=='..':
		    dir=os.path.split(dir)[0]
		else:
		    dir=os.path.join(dir,value)
	    self.configure(directory=dir)
	    self.fillit()

    def askfilename(self,directory=None,filter=None):
	"""The actual client function. Activates the dialog, and
	   returns only after a valid filename has been entered 
           (return value is that filename) or when canceled (return 
           value is None)"""
	if directory!=None:
	    self.configure(directory=directory)
	if filter!=None:
	    self.configure(filter=filter)
	self.fillit()
        self.canceled=1 # Needed for when user kills dialog window
	self.activate()
	if self.canceled:
	    return None
	else:
	    return self.component('filename').get()

    lastdir=""
    lastfilter=None
    lasttime=0
    def fillit(self):
	"""Get the directory list and show it in the two listboxes"""
        # Do not run unnecesarily
        if self.lastdir==self['directory'] and self.lastfilter==self['filter'] and self.lasttime>os.stat(self.lastdir)[8]:
            return
        self.lastdir=self['directory']
        self.lastfilter=self['filter']
        self.lasttime=time.time()
	dir=self['directory']
	if not dir:
	    dir=os.getcwd()
	dirs=['..']
	files=[]
        try:
            fl=os.listdir(dir)
            fl.sort()
        except os.error,arg:
            if arg[0] in (2,20):
                return
            raise
	for f in fl:
	    if os.path.isdir(os.path.join(dir,f)):
		dirs.append(f)
	    else:
		filter=self['filter']
		if not filter:
		    filter='*'
		if fnmatch.fnmatch(f,filter):
		    files.append(f)
	self.component('filenamebox').setlist(files)
	self.component('dirnamebox').setlist(dirs)
    
    def validate(self,filename):
	"""Validation function. Should return 1 if the filename is valid, 
           0 if invalid. May pop up dialogs to tell user why. Especially 
           suited to subclasses: i.e. only return 1 if the file does/doesn't 
           exist"""
	return 1

class FileDialogButtonClassFactory:
    def get(fn,filter='*'):
        """This returns a FileDialogButton class that will
        call the specified function with the resulting file.
        """
        class FileDialogButton(Tkinter.Button):
            # This is just an ordinary button with special colors.

            def __init__(self, master=None, cnf={}, **kw):
                '''when we get a file, we call fn(filename)'''
                self.fn = fn
                self.__toggle = 0
                apply(Tkinter.Button.__init__, (self, master, cnf), kw)
                self.configure(command=self.set)
            def set(self):
                fd = PmwFileDialog(self.master,filter=filter)
                fd.title('Please choose a file')
                n=fd.askfilename()
                if n is not None:
                    self.fn(n)
        return FileDialogButton
    get = staticmethod(get)

class DistMatTool:
    def setDistFilename(self,value):
        self.distfilename.setvalue(value)

    def setCaDistFilename(self,value):
        self.cadistfilename.setvalue(value)

    def __init__(self,app):
        parent = app.root
        self.parent = parent
        
        try:
                    self.min = stored.minValue
        except:
                        self.min = 0.0
            
        try:
                    self.mid = stored.midValue
        except:
                        self.mid = 0.0
            
        try:
                    self.max = stored.maxValue
        except:
                        self.max = 0.0
            
        self.incr = 0.01
        
        try:
                    self.distfilename = stored.distfilename
        except:
                        self.distfilename = 'MI.mat'

        try:
                    self.cadistfilename = stored.cadistfilename
        except:
                        self.cadistfilename = 'ca_dist.mat'

        string_distfilename = self.distfilename
        string_cadistfilename = self.cadistfilename

        self.cadistminvalue = 0.0
        try:
                    self.cadistmaxvalue = stored.cadistmaxvalue
        except:
                    self.cadistmaxvalue = 999.9

        self.cadistfilterFlag = False
        self.savestatFlag = False

        def quickFileValidation(s):
            if s == '': return Pmw.PARTIAL
            elif os.path.isfile(s): return Pmw.OK
            elif os.path.exists(s): return Pmw.PARTIAL
            else: return Pmw.PARTIAL
        
        self.dialog = Pmw.Dialog(parent,
                                 buttons = ('Show edges', 'Remove edges', 'Show neighbours', 'Update Cutoff Range', 'Exit showMIweb Tool'),
                                 title = 'PyMOL showMIweb Tool',
                                 command = self.execute)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        w = Tkinter.Label(self.dialog.interior(),
                                text = 'PyMOL showMIweb Tool\nAlessandro Pandini (C) 2009-14\n$Rev: 1531 $',
                                background = 'black',
                                foreground = 'white',
                                )
        w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)

        agroup = Pmw.Group(self.dialog.interior(),tag_text='Selection and filenames')
        agroup.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.object = Pmw.EntryField(agroup.interior(),
                                        labelpos='w',
                                        label_text='Object: ',
                                        value='all',
                                        )
        self.distfilename = Pmw.EntryField(agroup.interior(),
                                     labelpos='w',
                                     label_pyclass = FileDialogButtonClassFactory.get(self.setDistFilename),
                                     validate = {'validator':quickFileValidation,},
                                     value = string_distfilename,
                                     label_text = 'MI matrix filename:',
                                     )
        self.cadistfilename = Pmw.EntryField(agroup.interior(),
                                     labelpos='w',
                                     label_pyclass = FileDialogButtonClassFactory.get(self.setCaDistFilename),
                                     validate = {'validator':quickFileValidation,},
                                     value = string_cadistfilename,
                                     label_text = 'C-alpha distance matrix filename:',
                                     )
        bgroup = Pmw.Group(self.dialog.interior(),tag_text='Cutoffs')
        bgroup.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.cutoff = Pmw.Counter(bgroup.interior(),
                                        labelpos = 'w',
                                        label_text = 'MI Cutoff',
                                        label_justify = 'left',
                                        entryfield_value = round(self.mid * 100)/100,
                                        datatype = {'counter' : 'real', 'separator' : '.'},
                                        entryfield_validate = {'validator' : 'real',
                                                               'min' : self.min, 'max' : self.max,
                                                               'separator' : '.'},
                                        increment = self.incr)
        self.cadistmin = Pmw.Counter(bgroup.interior(),
                                        labelpos = 'w',
                                        label_text = 'min C-alpha distance',
                                        label_justify = 'left',
                                        entryfield_value = self.cadistminvalue,
                                        datatype = {'counter' : 'real', 'separator' : '.'},
                                        entryfield_validate = {'validator' : 'real',
                                                               'min' : self.cadistminvalue, 'max' : self.cadistmaxvalue,
                                                               'separator' : '.'},
                                        increment = self.incr * 10)
        self.cadistmax = Pmw.Counter(bgroup.interior(),
                                        labelpos = 'w',
                                        label_text = 'max C-alpha distance',
                                        label_justify = 'left',
                                        entryfield_value = round(self.cadistmaxvalue * 10 + 1)/10,
                                        datatype = {'counter' : 'real', 'separator' : '.'},
                                        entryfield_validate = {'validator' : 'real',
                                                               'min' : self.cadistminvalue, 'max' : self.cadistmaxvalue,
                                                               'separator' : '.'},
                                        increment = self.incr * 10)
        cgroup = Pmw.Group(self.dialog.interior(),tag_text='Representation')
        cgroup.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.reptype = Pmw.OptionMenu(cgroup.interior(),
                                      labelpos = 'w',
                                      label_text = 'Representation type',
                                      items = ['cartoon', 'ribbon', 'sticks', 'spheres'],
                                      )
        self.protpart = Pmw.OptionMenu(cgroup.interior(),
                                      labelpos = 'w',
                                      label_text = 'Protein part',
                                      items = ['all', 'backbone', 'C Alpha'],
                                      )
        self.colorscale = Pmw.OptionMenu(cgroup.interior(),
                                      labelpos = 'w',
                                      label_text = 'color scale',
                                      items = ['none', 'greyscale', 'colorscale'],
                                      )
        self.dashproportional = Pmw.OptionMenu(cgroup.interior(),
                                      labelpos = 'w',
                                      label_text = 'proportional dashes',
                                      items = ['no', 'yes'],
                                      )
        dgroup = Pmw.Group(self.dialog.interior(),tag_text='Neighbourhood analysis')
        dgroup.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.fragnum = Pmw.EntryField(dgroup.interior(),
                                        labelpos='w',
                                        label_text='Fragment number: ',
                                        value=1,
                                        )
        self.neighorder = Pmw.Counter(dgroup.interior(),
                                        labelpos = 'w',
                                        label_text = 'neighboroudhood order',
                                        label_justify = 'left',
                                        entryfield_value = 1,
                                        datatype = {'counter' : 'integer'},
                                        entryfield_validate = {'validator' : 'integer',
                                                               'min' : 1, 'max' : 3},
                                        increment = 1)
        egroup = Pmw.Group(self.dialog.interior(),tag_text='Options')
        egroup.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.savestatcheckbuttons = Pmw.RadioSelect(egroup.interior(),
                                      buttontype = 'checkbutton',
                                      labelpos = 'w',
                                      command = self.savestatcheckbuttoncallback,
                                      label_text = 'Save MI stats by distance'
                                      )   
        self.savestatcheckbuttons.pack(side = 'bottom', expand = 0, padx = 10, pady = 10) 
        self.savestatcheckbuttons.add("")
        #self.savestatcheckbuttons.invoke(Pmw.END)
        self.cadistcheckbuttons = Pmw.RadioSelect(egroup.interior(),
                                      buttontype = 'checkbutton',
                                      labelpos = 'w',
                                      command = self.cadistcheckbuttoncallback,
                                      label_text = 'Filter by distance'
                                      )   
        self.cadistcheckbuttons.pack(side = 'bottom', expand = 0, padx = 10, pady = 10) 
        self.cadistcheckbuttons.add("")
        #self.cadistcheckbuttons.invoke(Pmw.END)
        self.repcheckbuttons = Pmw.RadioSelect(egroup.interior(),
                                      buttontype = 'checkbutton',
                                      labelpos = 'w',
                                      command = self.repcheckbuttoncallback,
                                      label_text = 'Change representation'
                                      )   
        self.repcheckbuttons.pack(side = 'bottom', expand = 0, padx = 10, pady = 10) 
        self.repcheckbuttons.add("")
        self.repcheckbuttons.invoke(Pmw.END)

        for entry in (self.object,self.distfilename,self.cadistfilename,self.cutoff, self.cadistmin, self.cadistmax, self.reptype,self.protpart,self.colorscale,self.dashproportional, self.fragnum, self.neighorder):
            entry.pack(fill='x',padx=4,pady=1) 

        self.showAppModal()

    def repcheckbuttoncallback(self, tag, state):
        if state:
           self.changerepresentationFlag = True
        else:
           self.changerepresentationFlag = False

    def cadistcheckbuttoncallback(self, tag, state):
        if state:
           self.cadistfilterFlag = True
        else:
           self.cadistfilterFlag = False

    def savestatcheckbuttoncallback(self, tag, state):
        if state:
           self.savestatFlag = True
        else:
           self.savestatFlag = False

    def showAppModal(self):
        self.dialog.show()
        
    def changeRepresentation(self,name):
            cmd.bg_color('white') 
            cmd.hide('everything', name)
            cmd.color('green', name)
            if self.reptype.getvalue() == 'sticks':
                    part = ''
                    if self.protpart.getvalue() == 'backbone':
                        part = ' and name c+n+o+ca'
                    cmd.show(self.reptype.getvalue(), name + part)
            else:
                    if self.reptype.getvalue() == 'spheres':
                            part = ''
                            if self.protpart.getvalue() == 'backbone':
                                part = ' and name c+n+o+ca'
                            if self.protpart.getvalue() == 'C Alpha':
                                part = ' and name ca'
                            cmd.show(self.reptype.getvalue(), name + part)
                    else:
                            cmd.show(self.reptype.getvalue(), name)
        
    def read_distfile(self, filename):
            fin = open(filename.getvalue())
            lines = fin.readlines()
            fin.close()
            if len(lines) != len(lines[0].split()):
                print "Error: file %s is not a square matrix!" % filename
            else:
                mat = []
                nrow = len(lines)
                ncol = nrow
                for lidx in range(nrow):
                    values = [float(x) for x in lines[lidx].split()]
                    mat.append(values)
            return(mat)

    def create_FMN_colors(self):
        cmd.set_color("FMN001", [ 0.00,  0.00,  1.00])
        cmd.set_color("FMN002", [ 0.04,  0.02,  0.96])
        cmd.set_color("FMN003", [ 0.08,  0.05,  0.92])
        cmd.set_color("FMN004", [ 0.12,  0.08,  0.87])
        cmd.set_color("FMN005", [ 0.16,  0.10,  0.84])
        cmd.set_color("FMN006", [ 0.20,  0.13,  0.79])
        cmd.set_color("FMN007", [ 0.24,  0.16,  0.75])
        cmd.set_color("FMN008", [ 0.28,  0.18,  0.71])
        cmd.set_color("FMN009", [ 0.33,  0.21,  0.67])
        cmd.set_color("FMN010", [ 0.36,  0.24,  0.63])
        cmd.set_color("FMN011", [ 0.41,  0.26,  0.59])
        cmd.set_color("FMN012", [ 0.45,  0.29,  0.55])
        cmd.set_color("FMN013", [ 0.49,  0.31,  0.51])
        cmd.set_color("FMN014", [ 0.53,  0.34,  0.47])
        cmd.set_color("FMN015", [ 0.57,  0.37,  0.43])
        cmd.set_color("FMN016", [ 0.61,  0.40,  0.38])
        cmd.set_color("FMN017", [ 0.65,  0.42,  0.35])
        cmd.set_color("FMN018", [ 0.69,  0.45,  0.31])
        cmd.set_color("FMN019", [ 0.73,  0.47,  0.26])
        cmd.set_color("FMN020", [ 0.77,  0.50,  0.22])
        cmd.set_color("FMN021", [ 0.82,  0.53,  0.18])
        cmd.set_color("FMN022", [ 0.85,  0.55,  0.14])
        cmd.set_color("FMN023", [ 0.89,  0.58,  0.10])
        cmd.set_color("FMN024", [ 0.94,  0.60,  0.06])
        cmd.set_color("FMN025", [ 0.98,  0.63,  0.02])
        cmd.set_color("FMN026", [ 1.00,  0.63,  0.00])
        cmd.set_color("FMN027", [ 1.00,  0.60,  0.00])
        cmd.set_color("FMN028", [ 1.00,  0.58,  0.00])
        cmd.set_color("FMN029", [ 1.00,  0.55,  0.00])
        cmd.set_color("FMN030", [ 1.00,  0.53,  0.00])
        cmd.set_color("FMN031", [ 1.00,  0.50,  0.00])
        cmd.set_color("FMN032", [ 1.00,  0.47,  0.00])
        cmd.set_color("FMN033", [ 1.00,  0.45,  0.00])
        cmd.set_color("FMN034", [ 1.00,  0.42,  0.00])
        cmd.set_color("FMN035", [ 1.00,  0.40,  0.00])
        cmd.set_color("FMN036", [ 1.00,  0.37,  0.00])
        cmd.set_color("FMN037", [ 1.00,  0.34,  0.00])
        cmd.set_color("FMN038", [ 1.00,  0.31,  0.00])
        cmd.set_color("FMN039", [ 1.00,  0.29,  0.00])
        cmd.set_color("FMN040", [ 1.00,  0.26,  0.00])
        cmd.set_color("FMN041", [ 1.00,  0.24,  0.00])
        cmd.set_color("FMN042", [ 1.00,  0.21,  0.00])
        cmd.set_color("FMN043", [ 1.00,  0.18,  0.00])
        cmd.set_color("FMN044", [ 1.00,  0.16,  0.00])
        cmd.set_color("FMN045", [ 1.00,  0.13,  0.00])
        cmd.set_color("FMN046", [ 1.00,  0.10,  0.00])
        cmd.set_color("FMN047", [ 1.00,  0.08,  0.00])
        cmd.set_color("FMN048", [ 1.00,  0.05,  0.00])
        cmd.set_color("FMN049", [ 1.00,  0.02,  0.00])
        cmd.set_color("FMN050", [ 1.00,  0.00,  0.00])

    def execute(self, result):
        self.create_FMN_colors()
        if result == 'Show edges':
            try:
                cadistmat = self.read_distfile(self.cadistfilename)
                ncrow = len(cadistmat)
                nccol = ncrow
                self.cadistmaxvalue = 0
                camaxV = self.cadistmaxvalue
                stored.cadistmaxvalue = camaxV
                for i in range(ncrow):
                    for j in range(i,nccol):
                        if self.cadistmaxvalue < cadistmat[i][j]:
                            self.cadistmaxvalue = cadistmat[i][j]
                            camaxV = self.cadistmaxvalue
                            stored.cadistmaxvalue = camaxV
                print "Maximum C-alpha distance is %8.3f" % self.cadistmaxvalue
            except:
                print "C-alpha distance file %s not found." % self.cadistfilename.getvalue()
            mat = self.read_distfile(self.distfilename)
            nrow = len(mat)
            ncol = nrow
            cutoff =  float(self.cutoff.getvalue())
            try:
                if self.changerepresentationFlag:
                    self.changeRepresentation(self.object.getvalue())
            except:
                print "Error: cannot change representation of object %s!" % self.object.getvalue()
            try:
                for d in stored.distancelist:
                    cmd.delete(d)
            except:
                pass
            try:
                stored.distancelist = []
                distMItable = []
                for i in range(nrow):
                    for j in range(i,ncol):
                        distdisplayFlag = True
                        if self.cadistfilterFlag:
                            try:
                                if (cadistmat[i][j] > float(self.cadistmax.getvalue())) or (cadistmat[i][j] < float(self.cadistmin.getvalue())):
                                    distdisplayFlag = False
                            except:
                                pass
                        if (mat[i][j] > cutoff) and distdisplayFlag:
                            try:
                                distMItable.append((i, j, mat[i][j], cadistmat[i][j]))
                            except:
                                pass
                            cmd.distance('%s_%s' % (i + 1,j + 1), "name CA and resi %s" % (i + 1), "name CA and resi %s" % (j + 1))
                            cmd.hide('label', '%s_%s' % (i + 1,j + 1))
                            cmd.set('dash_gap', 0, '%s_%s' % (i + 1,j + 1))
                            cmd.set('dash_width', 4, '%s_%s' % (i + 1,j + 1))
                            if self.colorscale.getvalue() == 'greyscale':
                                cmd.color("grey%02d" % int((mat[i][j] - cutoff)/(self.max - cutoff) * 99), '%s_%s' % (i + 1,j + 1))
                            if self.colorscale.getvalue() == 'colorscale':
                                cmd.color("FMN%03d" % (int((mat[i][j] / self.max * 49)) + 1), '%s_%s' % (i + 1,j + 1))
                            if self.dashproportional.getvalue() == 'yes':
                                cmd.set('dash_radius', (mat[i][j] / self.max) * 0.25, '%s_%s' % (i + 1,j + 1))
                            stored.distancelist.append('%s_%s' % (i + 1,j + 1))
                if self.savestatFlag:
                    outfilename = "%s_I%4.2f_%04.1f-%04.1fang.dat" % (self.distfilename.getvalue()[:-4], cutoff, float(self.cadistmax.getvalue()), float(self.cadistmin.getvalue()))
                    fout = open(outfilename, 'w')
                    for d in distMItable:
                        fout.write("%5d\t%5d\t%f\t%f\n" % (d[0] + 1, d[1] + 1, d[2], d[3]))
                    fout.close()
            except:
                print "Error: inconsistency in the residue numbering!"
        else:
                if result == 'Update Cutoff Range':
                    try:
                        cadistmat = self.read_distfile(self.cadistfilename)
                        ncrow = len(cadistmat)
                        nccol = ncrow
                        self.cadistmaxvalue = 0
                        camaxV = self.cadistmaxvalue
                        stored.cadistmaxvalue = camaxV
                        for i in range(ncrow):
                            for j in range(i,nccol):
                                if self.cadistmaxvalue < cadistmat[i][j]:
                                    self.cadistmaxvalue = cadistmat[i][j]
                                    camaxV = self.cadistmaxvalue
                                    stored.cadistmaxvalue = camaxV
                        print "Maximum C-alpha distance is %8.3f" % self.cadistmaxvalue
                    except:
                        print "C-alpha distance file %s not found." % self.cadistfilename.getvalue()
                    mat = self.read_distfile(self.distfilename)
                    nrow = len(mat)
                    ncol = nrow
                    minV = self.min
                    maxV = self.max
                    for i in range(nrow):
                        for j in range(i,ncol):
                            if self.min > mat[i][j]:
                                stored.minValue = mat[i][j] #- self.incr
                                self.min = stored.minValue 
                                minV = mat[i][j]
                            if self.max < mat[i][j]:
                                stored.maxValue = mat[i][j] #+ self.incr
                                self.max = stored.maxValue
                                maxV = mat[i][j]
                    if minV >= 0.0 and self.min < 0.0:
                        self.min = 0.0
                    stored.midValue = self.max - (self.max - self.min) / 2
                    self.mid = stored.midValue
                    print "max MI: %8.3f" % maxV
                    print "min MI: %8.3f" % minV
                    stored.distfilename = self.distfilename.getvalue()
                    stored.cadistfilename = self.cadistfilename.getvalue()
                    self.dialog.withdraw()
                else:
                    if result == 'Remove edges':
                            try:
                                for d in stored.distancelist:
                                    cmd.delete(d)
                            except:
                                pass
                    else:
                        if result == 'Show neighbours':
                            try:
                                mat = self.read_distfile(self.distfilename)
                                nrow = len(mat)
                                ncol = nrow
                                if self.cadistfilterFlag:
                                    try:
                                        cadistmat = self.read_distfile(self.cadistfilename)
                                    except:
                                        print "C-alpha distance file %s not found." % self.cadistfilename.getvalue()
                                selectedFrag = int(self.fragnum.getvalue())
                                fragIndex = selectedFrag - 1
                                if selectedFrag > nrow:
                                    print "Warning: fragment index is greater than maximum value (%d)" % nrow
                                else:
                                    cmd.set("sphere_transparency", 0.5)
                                    cmd.select("frag_%d" % selectedFrag, "resi %s and name CA" % selectedFrag)
                                    cmd.show("spheres", "resi %d and name CA" % selectedFrag)
                                    cmd.color("white", "resi %d and name CA" % selectedFrag)
                                    reslistFirst = []
                                    resstringFirst = ""
                                    for k in range(ncol):
                                        distdisplayFlag = True
                                        if self.cadistfilterFlag:
                                            try:
                                                if (cadistmat[fragIndex][k] > float(self.cadistmax.getvalue())) or (cadistmat[fragIndex][k] < float(self.cadistmin.getvalue())):
                                                    distdisplayFlag = False
                                            except:
                                                pass
                                        if (mat[fragIndex][k] > float(self.cutoff.getvalue())) & distdisplayFlag:
                                            reslistFirst.append( (k + 1) )
                                            resstringFirst = resstringFirst + "%d+" % (k + 1)
                                    cmd.select("neigh_%d" % selectedFrag, "resi %s and name CA" % resstringFirst)
                                    cmd.show("spheres", "neigh_%d" % selectedFrag)
                                    cmd.color("red", "neigh_%d" % selectedFrag)
                                    if int(self.neighorder.getvalue()) > 1:
                                        reslistSecond = []
                                        resstringSecond = ""
                                        for resIdx in reslistFirst:  
                                            l = resIdx - 1
                                            for k in range(ncol):
                                                distdisplayFlag = True
                                                if self.cadistfilterFlag:
                                                    try:
                                                        if (cadistmat[l][k] > float(self.cadistmax.getvalue())) or (cadistmat[l][k] < float(self.cadistmin.getvalue())):
                                                            distdisplayFlag = False
                                                    except:
                                                        pass
                                                if (k + 1) in reslistFirst:
                                                    continue
                                                if k == fragIndex:
                                                    continue
                                                if (mat[l][k] > float(self.cutoff.getvalue())) & distdisplayFlag:
                                                    reslistSecond.append( (k + 1) )
                                                    resstringSecond = resstringSecond + "%d+" % (k + 1)
                                        cmd.select("neigh2_%d" % selectedFrag, "resi %s and name CA" % resstringSecond)
                                        cmd.show("spheres", "neigh2_%d" % selectedFrag)
                                        cmd.color("orange", "neigh2_%d" % selectedFrag)
                                    if int(self.neighorder.getvalue()) == 3:
                                        reslistThird = []
                                        resstringThird = ""
                                        for resIdx in (reslistFirst + reslistSecond):  
                                            l = resIdx - 1
                                            for k in range(ncol):
                                                if (k + 1) in (reslistFirst + reslistSecond):
                                                    continue
                                                if k == fragIndex:
                                                    continue
                                                distdisplayFlag = True
                                                if self.cadistfilterFlag:
                                                    try:
                                                        if (cadistmat[l][k] > float(self.cadistmax.getvalue())) or (cadistmat[l][k] < float(self.cadistmin.getvalue())):
                                                            distdisplayFlag = False
                                                    except:
                                                        pass
                                                if (mat[l][k] > float(self.cutoff.getvalue())) & distdisplayFlag:
                                                    reslistThird.append( (k + 1) )
                                                    resstringThird = resstringThird + "%d+" % (k + 1)
                                        if len(reslistThird) > 0:
                                            cmd.select("neigh3_%d" % selectedFrag, "resi %s and name CA" % resstringThird)
                                            cmd.show("spheres", "neigh3_%d" % selectedFrag)
                                            cmd.color("yellow", "neigh3_%d" % selectedFrag)
                                        else:
                                            print "No fragment identified as third neighbour"
                            except:
                                pass
                        else:
                                if __name__ == '__main__':
                                    self.parent.destroy()
                                else:
                                    self.dialog.withdraw()
