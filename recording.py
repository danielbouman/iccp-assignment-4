# -*- coding: utf-8 -*-
"""
Created on Tue May 12 12:16:47 2015

@author: Hanselman
"""

import matplotlib as mpl

class Recorder:

    def __init__(self, filename, *args):
        self.file = mpl.wave.open(filename,'wb')
        self.file.setnchannels(1)
        self.file.setsampwidth(128)
        self.file.setframerate(4*44e3)
        

    def write_sample(self,value):
        self.file.writeframes(value)

    def add_to_file(self,value):
        