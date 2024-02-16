from build.management.commands.base_build import Command as BaseBuild
# from django.db import connection
# from django.conf import settings
from django.db.models import F 


from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment, ProteinFamily, Species
# from structure.models import Structure, StructureModel, StructureComplexModel, StatsText, PdbData, StructureModelpLDDT, StructureType, StructureExtraProteins
# import structure.assign_generic_numbers_gpcr as as_gn
# from residue.models import Residue
# from common.models import WebResource, WebLink
# from common.definitions import G_PROTEIN_DISPLAY_NAME as g_prot_dict
# from signprot.models import SignprotComplex
# from contactnetwork.cube import compute_interactions
# from interaction.models import StructureLigandInteraction

from common.alignment import Alignment

from collections import OrderedDict
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
        gpcr_segments = ProteinSegment.objects.filter(proteinfamily='GPCR')
        cross_class_similarities = OrderedDict()
        proteins = None
        i = 1
        for gpcr_class in human_parent_gpcr_families[:-1]:
            del proteins
            cross_class_similarities2 = OrderedDict()
            # letting PostgreSQL handle the data storage via query cache, instead of storing in RAM with python.

            gpcr_class_proteins = Protein.objects.all().annotate(family_slug=F('family__slug')).filter(family_slug__startswith=gpcr_class.slug)
            # gpcr_class_proteins = gpcr_class_proteins[:10] #uncomment for test

            #gpcr_class_proteins = list(gpcr_class_proteins)
            for gpcr_class2 in human_parent_gpcr_families[i:]:
                gpcr_class2_proteins = Protein.objects.all().annotate(family_slug=F('family__slug')).filter(family_slug__startswith=gpcr_class2.slug)
                # gpcr_class2_proteins = gpcr_class2_proteins[:10] #uncomment for test
                proteins = gpcr_class_proteins | gpcr_class2_proteins
                proteins = proteins.order_by('family_slug')
                #proteins = gpcr_class_proteins + list(gpcr_class2_proteins)

                cs_alignment = Alignment()
                # cs_alignment.load_reference_protein(reference_protein)
                cs_alignment.load_proteins(proteins)
                cs_alignment.load_segments(gpcr_segments)
                cs_alignment.build_alignment()
                cs_alignment.remove_non_generic_numbers_from_alignment() 
                cs_alignment.calculate_similarity_matrix()

                # for key in cs_alignment.similarity_matrix.keys():
                #    print(key,cs_alignment.similarity_matrix[key])

                #print()

                entry_name_2_position = {}
                pos = 0
                for protein in cs_alignment.proteins:
                    entry_name_2_position[protein.protein.entry_name] = pos
                    pos += 1


                max_similarity = 0
                for protein1 in gpcr_class_proteins:
                    for protein2 in gpcr_class2_proteins:
                        pos = entry_name_2_position[protein2.entry_name]
                        similarity_value = int(cs_alignment.similarity_matrix[protein1.entry_name]['values'][pos][0])
                        if  similarity_value > max_similarity:
                            max_similarity = similarity_value
                        #print(protein1,protein2,str(similarity_value))
                cross_class_similarities2[gpcr_class2.name] = max_similarity
            
            cross_class_similarities[gpcr_class.name] = cross_class_similarities2
            i += 1


        # for key in  cross_class_similarities:

        #     print(key,cross_class_similarities[key])

        matrix = OrderedDict()
        human_parent_gpcr_families_names = [gpcr_class.name for gpcr_class in human_parent_gpcr_families]
        for key in human_parent_gpcr_families_names:
            row = []
            for key2 in human_parent_gpcr_families_names:
                if key == key2:
                    row.append('-')
                    continue
                if key in cross_class_similarities:
                    if key2 in cross_class_similarities[key]:
                        row.append(str(cross_class_similarities[key][key2]))
                    else:
                        row.append(str(cross_class_similarities[key2][key]))
                else:
                    row.append(str(cross_class_similarities[key2][key]))
            matrix[key] = row

        headers = matrix.keys()
        print(', '.join(list(headers)))
        for key in  matrix:

            print(key,matrix[key])
             



