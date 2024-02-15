from build.management.commands.base_build import Command as BaseBuild
# from django.db import connection
# from django.conf import settings
from django.db.models import F 


from protein.models import Protein, ProteinFamily, Species, ProteinState, ProteinConformation
# from structure.models import Structure, StructureModel, StructureComplexModel, StatsText, PdbData, StructureModelpLDDT, StructureType, StructureExtraProteins
# import structure.assign_generic_numbers_gpcr as as_gn
# from residue.models import Residue
# from common.models import WebResource, WebLink
# from common.definitions import G_PROTEIN_DISPLAY_NAME as g_prot_dict
# from signprot.models import SignprotComplex
# from contactnetwork.cube import compute_interactions
# from interaction.models import StructureLigandInteraction

from common.alignment import Alignment

import Bio.PDB as PDB
import os
import logging
# import zipfile
# import shutil
from datetime import datetime, date
# import time


starttime = datetime.now()
logger = logging.getLogger('class_similarity')
hdlr = logging.FileHandler('./logs/class_similarity.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.INFO)
structure_path = './structure/'
pir_path = os.sep.join([structure_path, 'PIR'])

classless_parent_gpcr_slugs = set(('008')) # name = Other GPCRs

build_date = date.today()

import warnings
warnings.filterwarnings("ignore")

class Command(BaseBuild):
    help = 'Build cross-class similarity matrices and stores them in GPCRdb database.'

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        # parser.add_argument('-f', help='Specify file name to be uploaded to GPCRdb', default=False, type=str, nargs='+')
        # parser.add_argument('-c', help='Upload only complex models to GPCRdb', default=False, action='store_true')
        # parser.add_argument('--purge', help='Purge existing entries in GPCRdb', default=False, action='store_true')
        # parser.add_argument('--copy_from_data', help='Copy and rename AF models from data to upload folder', default=False, action='store_true')
                
    def get_parent_gpcr_families(self):
        parent_family = ProteinFamily.objects.get(slug='000') 
        parent_gpcr_families = ProteinFamily.objects.filter(parent_id=parent_family.pk, slug__startswith='0').exclude(pk=parent_family.pk)
        return sorted(parent_gpcr_families,key=lambda f: (int(f.slug)))
    
    def filter_out_non_human_parent_gpcr_families(self,parent_gpcr_families):
        parent_gpcr_families_ids = [f.id for f in parent_gpcr_families]
        
        q = Protein.objects.filter(family__in=parent_gpcr_families,species=Species.objects.get(common_name__iexact='Human'))
        human_parent_gpcr_families_ids = q.annotate(family_pk=F('family__pk')).distinct('family_pk').values_list('family_pk',flat=True)
        return [f for f in parent_gpcr_families if f.id in human_parent_gpcr_families_ids]  

    def handle(self, *args, **options):
        
        parent_families = self.get_parent_gpcr_families()
        human_parent_gpcr_families = self.filter_out_non_human_parent_gpcr_families(parent_families)
        for gpcr_class in human_parent_gpcr_families:
             print(gpcr_class)
             
        # cs_alignment = Alignment()
        # cs_alignment.load_reference_protein(reference_protein)
        # cs_alignment.load_proteins([i.protein_conformation.protein for i in list(Homology_model.similarity_table.keys())])


