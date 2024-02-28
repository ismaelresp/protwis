from build.management.commands.base_build import Command as BaseBuild
# from django.db import connection
# from django.conf import settings
from django.db.models import F 


from protein.models import Protein, ProteinConformation, ProteinAnomaly, ProteinState, ProteinSegment, ProteinFamily, Species
from alignment.models import ClassSimilarity
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
import sys
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

matrix_header_max_length = 20

build_date = date.today()

import warnings
warnings.filterwarnings("ignore")

class Command(BaseBuild):
    help = 'Build cross-class similarity matrices and stores them into GPCRdb database purging first the existing similirity entries.'
    

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose', help='Prints progress in stdout.', default=False, action='store_true')
        parser.add_argument('--print_results', help='Prints results in stdout.', default=False, action='store_true')
       
    def get_parent_gpcr_families(self):
        parent_family = ProteinFamily.objects.get(slug='000') 
        parent_gpcr_families = ProteinFamily.objects.filter(parent_id=parent_family.pk, slug__startswith='0').exclude(pk=parent_family.pk)
        return sorted(parent_gpcr_families,key=lambda f: (int(f.slug)))
    
    def get_human_species(self):
        return Species.objects.get(common_name__iexact='Human')
    
    def get_yeast_species(self):
        return Protein.objects.filter(entry_name__iendswith='_yeast')[0].species
    
    def filter_out_non_species_parent_gpcr_families(self,parent_gpcr_families,species):
        new_parent_gpcr_families_slugs_set = set()
        species2 = species
        try:
            species_iterator = iter(species)
        except TypeError as te:
             species2 = [species]
        parent_gpcr_families_slugs = []
        for family in parent_gpcr_families:
            parent_gpcr_families_slugs.append(family.slug)

        for slug in parent_gpcr_families_slugs:
            q = Protein.objects.annotate(family_slug=F('family__slug')).filter(family_slug__startswith=slug,species__in=species2)
            if q.exists():
                new_parent_gpcr_families_slugs_set.add(slug)
        return [f for f in parent_gpcr_families if f.slug in new_parent_gpcr_families_slugs_set]
        
    
    def filter_out_non_human_parent_gpcr_families(self,parent_gpcr_families):
        return self.filter_out_non_species_parent_gpcr_families(parent_gpcr_families,self.get_human_species())

    def filter_out_non_yeast_parent_gpcr_families(self,parent_gpcr_families):
        return self.filter_out_non_species_parent_gpcr_families(parent_gpcr_families,self.get_yeast_species())
    
    def handle(self, *args, **options):
        initial_step1 = 380 #If alignment fails, please, set this to a lower value
        initial_step2 = 380 #If alignment fails, please, set this to a lower value


        parent_families = self.get_parent_gpcr_families()
        human_parent_gpcr_families = self.filter_out_non_human_parent_gpcr_families(parent_families)
        human_species = self.get_human_species()
        yeast_parent_gpcr_families = self.filter_out_non_yeast_parent_gpcr_families(parent_families)
        yeast_non_human_parent_gpcr_families_set = set(yeast_parent_gpcr_families) - set(human_parent_gpcr_families)
        yeast_non_human_parent_gpcr_families = [gpcr_family for gpcr_family in yeast_parent_gpcr_families if gpcr_family in yeast_non_human_parent_gpcr_families_set]
        yeast_species = self.get_yeast_species()
        gpcr_segments = ProteinSegment.objects.filter(proteinfamily='GPCR')

        cross_class_most_similar_gpcr_pairs = {}
        cross_class_similarities = OrderedDict()
        proteins = None
        i = 1



        step1=int(initial_step1) #If alignment fails, please, set this to a lower value
        step2=int(initial_step2) #If alignment fails, please, set this to a lower value
        step_halved = False 


        human_parent_gpcr_families_protein = {}
        human_parent_gpcr_families_protein_num = {}
        for gpcr_class in human_parent_gpcr_families:
            gpcr_class_proteins = Protein.objects.all().annotate(family_slug=F('family__slug')).filter(species=human_species,family_slug__startswith=gpcr_class.slug)
            gpcr_class_proteins = gpcr_class_proteins.exclude(accession=None).order_by('family_slug','entry_name')
            human_parent_gpcr_families_protein[gpcr_class] = list(gpcr_class_proteins)
            human_parent_gpcr_families_protein_num[gpcr_class] = len(gpcr_class_proteins)

        yeast_non_human_parent_gpcr_families_protein = {}
        yeast_non_human_parent_gpcr_families_protein_num = {} 
        for gpcr_class in yeast_non_human_parent_gpcr_families:
            gpcr_class_proteins = Protein.objects.all().annotate(family_slug=F('family__slug')).filter(species=yeast_species,family_slug__startswith=gpcr_class.slug)
            gpcr_class_proteins = gpcr_class_proteins.exclude(accession=None).order_by('family_slug','entry_name')
            yeast_non_human_parent_gpcr_families_protein[gpcr_class] = list(gpcr_class_proteins)
            yeast_non_human_parent_gpcr_families_protein_num[gpcr_class] = len(gpcr_class_proteins)

        selected_parent_gpcr_families = human_parent_gpcr_families + yeast_non_human_parent_gpcr_families
        selected_parent_gpcr_families_protein = {}
        selected_parent_gpcr_families_protein_num = {}
        for gpcr_class in selected_parent_gpcr_families:
            if gpcr_class in human_parent_gpcr_families_protein and gpcr_class in yeast_non_human_parent_gpcr_families_protein:
                selected_parent_gpcr_families_protein[gpcr_class] = human_parent_gpcr_families_protein[gpcr_class] + yeast_non_human_parent_gpcr_families_protein[gpcr_class]
                selected_parent_gpcr_families_protein_num[gpcr_class] = human_parent_gpcr_families_protein_num[gpcr_class] + yeast_non_human_parent_gpcr_families_protein_num[gpcr_class]
            elif gpcr_class in human_parent_gpcr_families_protein:
                selected_parent_gpcr_families_protein[gpcr_class] = human_parent_gpcr_families_protein[gpcr_class]
                selected_parent_gpcr_families_protein_num[gpcr_class] = human_parent_gpcr_families_protein_num[gpcr_class]
            elif gpcr_class in yeast_non_human_parent_gpcr_families_protein:
                selected_parent_gpcr_families_protein[gpcr_class] = yeast_non_human_parent_gpcr_families_protein[gpcr_class]
                selected_parent_gpcr_families_protein_num[gpcr_class] = yeast_non_human_parent_gpcr_families_protein_num[gpcr_class]


        for gpcr_class in selected_parent_gpcr_families[:-1]:
            while_loop_continue = False
            while True:
                cross_class_similarities2 = OrderedDict()
                cross_class_most_similar_gpcr_pairs2 = {}
                for clim in range(0,selected_parent_gpcr_families_protein_num[gpcr_class],step1):
                    
                    gpcr_class_proteins = selected_parent_gpcr_families_protein[gpcr_class]
                    gpcr_class_proteins = gpcr_class_proteins[clim:clim+step1]
                    # gpcr_class_proteins = gpcr_class_proteins[:10] #uncomment for testing
                    
                    for gpcr_class2 in selected_parent_gpcr_families[i:]:
                        
                        for clim2 in range(0,selected_parent_gpcr_families_protein_num[gpcr_class2],step2):
                            if options['verbose']:
                                print(gpcr_class,"from:"+str(clim+1),"to:"+str(clim+step1),"(of:{})".format(selected_parent_gpcr_families_protein_num[gpcr_class]),\
                                    'vs',gpcr_class2,"from:"+str(clim2+1),"to:"+str(clim2+step2),"(of:{})".format(selected_parent_gpcr_families_protein_num[gpcr_class2]))
                            gpcr_class2_proteins = selected_parent_gpcr_families_protein[gpcr_class2]
                            gpcr_class2_proteins = gpcr_class2_proteins[clim2:clim2+step2]
                            # gpcr_class2_proteins = gpcr_class2_proteins[:10] #uncomment for testing
                            proteins = gpcr_class_proteins + gpcr_class2_proteins

                            cs_alignment = Alignment()
                            cs_alignment.load_proteins(proteins)
                            cs_alignment.load_segments(gpcr_segments)
                            build_alignment_return_value = cs_alignment.build_alignment()
                            if build_alignment_return_value == "Too large":
                                print('Alignment too large. Retrying...', file=sys.stderr)
                                while_loop_continue = True
                                break
                            cs_alignment.remove_non_generic_numbers_from_alignment() 
                            cs_alignment.calculate_similarity_matrix()

                            entry_name_2_position = {}
                            pos = 0
                            for protein in cs_alignment.proteins:
                                entry_name_2_position[protein.protein.entry_name] = pos
                                pos += 1


                            max_similarity = 0
                            for protein1 in gpcr_class_proteins:
                                for protein2 in gpcr_class2_proteins:
                                    pos = entry_name_2_position[protein1.entry_name]
                                    similarity_value = int(cs_alignment.similarity_matrix[protein2.entry_name]['values'][pos][0])
                                    if  similarity_value > max_similarity:
                                        max_similarity = similarity_value
                                        protein_tuple = (protein1,protein2)
                                        
                            if gpcr_class2.name in cross_class_similarities2:
                                if cross_class_similarities2[gpcr_class2.name] < max_similarity:
                                    cross_class_similarities2[gpcr_class2.name] = max_similarity
                                    cross_class_most_similar_gpcr_pairs2[gpcr_class2] = protein_tuple
                            else:
                                cross_class_similarities2[gpcr_class2.name] = max_similarity
                                cross_class_most_similar_gpcr_pairs2[gpcr_class2] = protein_tuple
                        if while_loop_continue:
                            break
                    del proteins
                    del gpcr_class_proteins
                    if while_loop_continue:
                        break
                    selected_parent_gpcr_families_protein[gpcr_class][clim:clim+step1] = [None for i1 in range(clim,clim+step1)]
                if while_loop_continue:
                    if options['verbose']: print("Halving step1 and step2...") 
                    step1 = step1 // 2
                    step2 = step2 // 2
                    step_halved = True
                    while_loop_continue = False
                    if options['verbose']: print("Retrying last class similarity computation with the new steps...")
                    continue
                elif step_halved:
                    print("Restoring initial step1 and step2...")
                    step1=initial_step1
                    step2=initial_step2
                    step_halved = False
                cross_class_similarities[gpcr_class.name] = cross_class_similarities2
                cross_class_most_similar_gpcr_pairs[gpcr_class] = cross_class_most_similar_gpcr_pairs2
                del selected_parent_gpcr_families_protein[gpcr_class]
                i += 1
                break    

        # for key in  cross_class_similarities:
        #   print(key,cross_class_similarities[key])
        if options['print_results']:
            print("[Protein_pairs]")
            for gpcr_class in cross_class_most_similar_gpcr_pairs:      
                for gpcr_class2 in cross_class_most_similar_gpcr_pairs[gpcr_class]:
                    protein1,protein2 = cross_class_most_similar_gpcr_pairs[gpcr_class][gpcr_class2]
                    print(';'.join((gpcr_class.name,gpcr_class2.name,protein1.entry_name,protein2.entry_name)))
            print("[End Protein_pairs]")

            cross_class_similarity_matrix = OrderedDict()
            selected_parent_gpcr_families_names = [gpcr_class.name for gpcr_class in selected_parent_gpcr_families]
            for key in selected_parent_gpcr_families_names:
                row = []
                for key2 in selected_parent_gpcr_families_names:
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
                cross_class_similarity_matrix[key] = row

            print("[Matrix]")
            headers = cross_class_similarity_matrix.keys()
            print(';'+';'.join(list(headers)))
            for key in  cross_class_similarity_matrix:
                print(("{};"+';'.join(["{}" for i in range(0,len(cross_class_similarity_matrix[key]))])).format(key,*cross_class_similarity_matrix[key]))
            print("[End Matrix]")
            if options['verbose']: 
                print("[Human readable Protein_pairs]")
                gpcr_class_name_max_len = max([len(g.name) for g in cross_class_most_similar_gpcr_pairs])
                for gpcr_class in cross_class_most_similar_gpcr_pairs:      
                    for gpcr_class2 in cross_class_most_similar_gpcr_pairs[gpcr_class]:
                        protein1,protein2 = cross_class_most_similar_gpcr_pairs[gpcr_class][gpcr_class2]
                        print('{:{w}s} vs {:{w}s}'.format(gpcr_class.name[:matrix_header_max_length],gpcr_class2.name[:matrix_header_max_length],\
                                                        w=matrix_header_max_length),'{:>{w}s} vs {:>{w}s}'.format(protein1.entry_name,protein2.entry_name,\
                                                                                                                w=25))
                print("[End Human readable Protein_pairs]")

                print("[Human readable Matrix]")
                headers = cross_class_similarity_matrix.keys()
                print('; '+'; '.join(list(headers)))
                for key in  cross_class_similarity_matrix:
                    nkey = key[:matrix_header_max_length]
                    print(("{:{header_max_length}s} "+' '.join(["{: >3}" for i in range(0,len(cross_class_similarity_matrix[key]))]))\
                        .format(nkey,*cross_class_similarity_matrix[key],header_max_length=matrix_header_max_length))
                print("[End Human readable Matrix]")

        if options['verbose']: print('Truncating old table...')
        ClassSimilarity.custom_objects.truncate_table() #custom_objects is a custom manager, see ClassSimilarity model definition
        if options['verbose']: print('Saving new table...')
        class_similarity_list = []
        for gpcr_class in cross_class_most_similar_gpcr_pairs:      
            for gpcr_class2 in cross_class_most_similar_gpcr_pairs[gpcr_class]:
                protein1,protein2 = cross_class_most_similar_gpcr_pairs[gpcr_class][gpcr_class2]
                class_similarity_list.append(ClassSimilarity(protein_family1=gpcr_class,protein_family2=gpcr_class2,protein1=protein1,protein2=protein2,similarity=cross_class_similarities[gpcr_class.name][gpcr_class2.name]))
        ClassSimilarity.objects.bulk_create(class_similarity_list)
        if options['verbose']: print('Done.')

     





             



