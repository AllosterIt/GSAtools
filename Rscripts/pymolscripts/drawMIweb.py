#!/usr/bin/env python
'''
$Id: drawMIweb.py 1390 2013-04-16 00:42:22Z apandini $

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
import Tkinter
from Tkinter import *
import Pmw
from pymol import cmd
from pymol import stored

def __init__(self):
    self.menuBar.addmenuitem(   'Plugin', 'command',
                                'drawMIweb',
                                label = 'drawMIweb',
                                command = lambda s=self : DistMatTool(s))

class DistMatTool:
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
            
        self.incr = 0.1
        
        try:
                    self.filename = stored.filename
        except:
                        self.filename = 'MI.mat'
        
        self.dialog = Pmw.Dialog(parent,
                                 buttons = ('Show edges', 'Remove edges', 'Update Cutoff Range', 'Exit drawMIweb Tool'),
                                 title = 'PyMOL drawMIweb Tool',
                                 command = self.execute)
        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        w = Tkinter.Label(self.dialog.interior(),
                                text = 'PyMOL drawMIweb Tool\nAlessandro Pandini (C) 2009-13',
                                background = 'black',
                                foreground = 'white',
                                )
        w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)

        group = Pmw.Group(self.dialog.interior(),tag_text='Main options')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        self.object = Pmw.EntryField(group.interior(),
                                        labelpos='w',
                                        label_text='Object: ',
                                        value='all',
                                        )
        self.distfilename = Pmw.EntryField(group.interior(),
                                        labelpos='w',
                                        label_text='MI matrix filename',
                                        value= self.filename,
                                        )
        self.cutoff = Pmw.Counter(group.interior(),
                                        labelpos = 'w',
                                        label_text = 'MI Cutoff',
                                        label_justify = 'left',
                                        entryfield_value = self.mid,
                                        datatype = {'counter' : 'real', 'separator' : '.'},
                                        entryfield_validate = {'validator' : 'real',
                                                               'min' : self.min, 'max' : self.max,
                                                               'separator' : '.'},
                                        increment = self.incr)
        self.reptype = Pmw.OptionMenu(group.interior(),
                                      labelpos = 'w',
                                      label_text = 'Representation type',
                                      items = ['cartoon', 'ribbon', 'sticks', 'spheres'],
                                      )
        self.protpart = Pmw.OptionMenu(group.interior(),
                                      labelpos = 'w',
                                      label_text = 'Protein part',
                                      items = ['all', 'backbone', 'C Alpha'],
                                      )
        self.colorscale = Pmw.OptionMenu(group.interior(),
                                      labelpos = 'w',
                                      label_text = 'color scale',
                                      items = ['none', 'greyscale', 'colorscale'],
                                      )
        self.dashproportional = Pmw.OptionMenu(group.interior(),
                                      labelpos = 'w',
                                      label_text = 'proportional dashes',
                                      items = ['no', 'yes'],
                                      )

        for entry in (self.object,self.distfilename,self.cutoff,self.reptype,self.protpart,self.colorscale,self.dashproportional):
            entry.pack(fill='x',padx=4,pady=1) 

        self.showAppModal()

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
            mat = self.read_distfile(self.distfilename)
            nrow = len(mat)
            ncol = nrow
            cutoff =  float(self.cutoff.getvalue())
            try:
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
                for i in range(nrow):
                    for j in range(i,ncol):
                        if mat[i][j] > cutoff:
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
            except:
                print "Error: inconsistency in the residue numbering!"
        else:
                if result == 'Update Cutoff Range':
                    mat = self.read_distfile(self.distfilename)
                    nrow = len(mat)
                    ncol = nrow
                    minV = self.min
                    maxV = self.max
                    for i in range(nrow):
                        for j in range(i,ncol):
                            if self.min > mat[i][j]:
                                stored.minValue = mat[i][j] - self.incr
                                self.min = stored.minValue 
                                minV = mat[i][j]
                            if self.max < mat[i][j]:
                                stored.maxValue = mat[i][j] + self.incr
                                self.max = stored.maxValue
                                maxV = mat[i][j]
                    if minV >= 0.0 and self.min < 0.0:
                        self.min = 0.0
                    stored.midValue = self.max - (self.max - self.min) / 2
                    self.mid = stored.midValue
                    print "max %8.3f" % maxV
                    print "min %8.3f" % minV
                    stored.filename = self.distfilename.getvalue()
                    self.dialog.withdraw()
                else:
                    if result == 'Remove edges':
                            try:
                                for d in stored.distancelist:
                                    cmd.delete(d)
                            except:
                                pass
                    else:
                            if __name__ == '__main__':
                                self.parent.destroy()
                            else:
                                self.dialog.withdraw()
