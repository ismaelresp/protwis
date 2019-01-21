from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command

from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from django.db import IntegrityError
from contactnetwork.cube import compute_interactions
from contactnetwork.models import *
import contactnetwork.interaction as ci

from residue.models import ResidueGenericNumber, ResidueNumberingScheme, Residue, ResidueGenericNumberEquivalent


from structure.models import (Structure, StructureType, StructureSegment, StructureStabilizingAgent,PdbData,
    Rotamer, StructureSegmentModeling, StructureCoordinates, StructureCoordinatesDescription, StructureEngineering,
    StructureEngineeringDescription, Fragment)


import os
import yaml
from interaction.views import runcalculation,parsecalculation
from multiprocessing import Queue, Process, Value, Lock

class Command(BaseCommand):

    help = "Output all uniprot mappings"

    def prepare_input(self, proc, items, iteration=1):
        q = Queue()
        procs = list()
        num_items = len(items)
        num = Value('i', 0)
        lock = Lock()

        if not num_items:
            return False

        # make sure not to use more jobs than proteins (chunk size will be 0, which is not good)
        if proc > num_items:
            proc = num_items

        chunk_size = int(num_items / proc)
        connection.close()
        for i in range(0, proc):
            first = chunk_size * i
            if i == proc - 1:
                last = False
            else:
                last = chunk_size * (i + 1)
    
            p = Process(target=self.main_func, args=([(first, last), iteration,num,lock]))
            procs.append(p)
            p.start()

        for p in procs:
            p.join()

    def purge_contact_network(self,s):

        ii = InteractingResiduePair.objects.filter(
            referenced_structure=s
        ).all().delete()

    def build_contact_network(self,s,pdb_code):
        interacting_pairs, distances  = compute_interactions(pdb_code, save_to_db=True)


    def handle(self, *args, **options):

        self.ss = Structure.objects.filter(refined=False).all()
        self.structure_data_dir = os.sep.join([settings.DATA_DIR, 'structure_data', 'structures'])
        self.prepare_input(4, self.ss)

        # for s in Structure.objects.filter(refined=False).all():
        #   self.purge_contact_network(s)
        #   self.build_contact_network(s,s.pdb_code.index)

    def main_func(self, positions, iteration,count,lock):
        # filenames
        # if not positions[1]:
        #     filenames = self.filenames[positions[0]:]
        # else:
        #     filenames = self.filenames[positions[0]:positions[1]]
        ss = self.ss
        while count.value<len(ss):
            with lock:
                if count.value<len(ss):
                    s = ss[count.value]
                    count.value +=1
                else:
                    break 

            source_file_path = os.sep.join([self.structure_data_dir, s.pdb_code.index.upper() + ".yaml"])
            if os.path.isfile(source_file_path):
                with open(source_file_path, 'r') as f:
                    sd = yaml.load(f)
                    
            peptide_chain = ""
            if 'ligand' in sd and sd['ligand'] and sd['ligand']!='None':
                if isinstance(sd['ligand'], list):
                    ligands = sd['ligand']
                else:
                    ligands = [sd['ligand']]
                for ligand in ligands:
                    peptide_chain = ""
                    if 'chain' in ligand:
                        peptide_chain = ligand['chain']
            
            # self.purge_contact_network(s)
            self.build_contact_network(s,s.pdb_code.index)
            print(s,"Contact Network",time.time()-current)
            current = time.time()
            runcalculation(s.pdb_code.index,peptide_chain)
            parsecalculation(s.pdb_code.index,False)
            print(s,"Ligand Interactions",time.time()-current)