# -*- coding: utf-8 -*-
class Material():
    def __init__(self,elastic,stress,pla):
        self.elastic=elastic
        self.pla=bool(pla)
        if pla:
            self.stress=stress
class Section():
    def __init__(self,mname,elname):
        self.material_name=mname
        self.elset_name=elname
class Object():
    def __init__(self):
        self.nset={}
        self.elset={}
        self.materials={}
        self.sections={}
        self.ngen={}
        self.elgen={}
        self.sections={}
    def add_nset(self,nname,nodes,ngen):
        self.nset[nname]=nodes
        self.ngen[nname]=bool(ngen)
    def add_elset(self,elname,elements,elgen):
        self.elset[elname]=elements
        self.elgen[elname]=bool(elgen)
    def del_elset(self,elname):
        self.elset.pop(elname)
    def add_material(self,mname,elastic,stress,pla):
        self.materials[mname]=Material(elastic,stress,pla)
    def change_material(self,mname,elastic):
        self.materials[mname].elastic=elastic
    def add_section(self,secname,mname,elname):
        self.sections[secname]=Section(mname, elname)
        
        