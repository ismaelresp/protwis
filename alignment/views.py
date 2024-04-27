from django.shortcuts import render, redirect
from django.conf import settings
from django.http import HttpResponse
from django.views.generic import TemplateView
from django.db.models import Case, When
from django.core.cache import cache
from django.core.cache import caches
from django.utils.html import escape
from django.core.files.base import File

try:
    cache_alignment = caches['alignments']
except:
    cache_alignment = cache

from alignment.functions import get_proteins_from_selection
from common import definitions
from common.selection import Selection, SelectionItem
from common.views import AbsTargetSelection, AbsTargetSelectionTable
from common.views import AbsSegmentSelection
from common.views import AbsMiscSelection
from structure.functions import BlastSearch
from protwis.context_processors import site_title

# from common.alignment_SITE_NAME import Alignment
Alignment = getattr(__import__('common.alignment_' + settings.SITE_NAME, fromlist=['Alignment']), 'Alignment')
from protein.models import Protein, ProteinSegment, ProteinFamily, ProteinSet
from residue.models import ResidueNumberingScheme, ResiduePositionSet
from alignment.models import ClassSimilarity, ClassSimilarityTie, ClassSimilarityType, ClassRepresentativeSpecies
from seqsign.sequence_signature import SequenceSignature, signature_score_excel
from common.definitions import CLASSLESS_PARENT_GPCR_SLUGS

from collections import OrderedDict
from copy import deepcopy
import hashlib
import inspect
from io import BytesIO
import itertools
import json
import numpy as np
import os
import xlsxwriter
from xlsxwriter.utility import xl_range_abs
import xlrd
import re

strain_re = re.compile(r'\bstrain\b', flags=re.I)
class_fungal_re = re.compile(r'(\(Ste2-like)(\s+)(fungal)(\s+)(pheromone\))', flags=re.I)
class_fullname_re = re.compile(r'(Class\s+\w+)(\s+)(\(.*\))', flags=re.I)

# class TargetSelection(AbsTargetSelection):
#     step = 1
#     number_of_steps = 2
#     filter_tableselect = False
#     docs = 'sequences.html#structure-based-alignments'
#     selection_boxes = OrderedDict([
#         ('reference', False),
#         ('targets', True),
#         ('segments', False),
#     ])
#     buttons = {
#         'continue': {
#             'label': 'Continue to next step',
#             'url': '/alignment/segmentselection',
#             'color': 'success',
#         },
#     }

class TargetSelection(AbsTargetSelectionTable):
    step = 1
    number_of_steps = 2
    filter_tableselect = False
    docs = 'sequences.html#structure-based-alignments'
    title = "SELECT RECEPTORS"
    description = 'Select receptors in the table (below) or browse the classification tree (right). You can select entire' \
        + ' families or individual receptors.\n\nOnce you have selected all your receptors, click the green button.'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Next',
            'onclick': "submitSelection('/alignment/segmentselection');",
            'color': 'success',
        },
    }

class TargetSelectionGprotein(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    psets = False
    filters = True
    filter_gprotein = True

    docs = 'sequences.html#structure-based-alignments'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselectiongprot',
            'color': 'success',
        },
    }
    try:
        if ProteinFamily.objects.filter(slug="100_001").exists():
            ppf = ProteinFamily.objects.get(slug="100_001")
            pfs = ProteinFamily.objects.filter(parent=ppf.id)
            ps = Protein.objects.filter(family=ppf)

            tree_indent_level = []
            action = 'expand'
            # remove the parent family (for all other families than the root of the tree, the parent should be shown)
            del ppf
    except Exception as e:
        pass

class TargetSelectionArrestin(AbsTargetSelection):
    step = 1
    number_of_steps = 2
    psets = False
    filters = True
    filter_gprotein = True

    docs = 'sequences.html#structure-based-alignments'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', False),
    ])
    buttons = {
        'continue': {
            'label': 'Continue to next step',
            'url': '/alignment/segmentselectionarrestin',
            'color': 'success',
        },
    }
    try:
        if ProteinFamily.objects.filter(slug="200_000").exists():
            ppf = ProteinFamily.objects.get(slug="200_000")
            pfs = ProteinFamily.objects.filter(parent=ppf.id)
            ps = Protein.objects.filter(family=ppf)

            tree_indent_level = []
            action = 'expand'
            # remove the parent family (for all other families than the root of the tree, the parent should be shown)
            del ppf
    except Exception as e:
        pass

class SegmentSelection(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }

class SegmentSelectionGprotein(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    description = 'Select sequence segments in the middle column for G proteins. You can expand every structural element and select individual' \
        + ' residues by clicking on the down arrows next to each helix, sheet or loop.\n\n You can select the full sequence or show all structured regions at the same time.\n\nSelected segments will appear in the' \
        + ' right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green' \
        + ' button.'

    template_name = 'common/segmentselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }

    position_type = 'gprotein'
    rsets = ResiduePositionSet.objects.filter(name__in=['Gprotein Barcode', 'YM binding site']).prefetch_related('residue_position')

    ss = ProteinSegment.objects.filter(partial=False, proteinfamily='Alpha').prefetch_related('generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')


class SegmentSelectionArrestin(AbsSegmentSelection):
    step = 2
    number_of_steps = 2
    docs = 'sequences.html#structure-based-alignments'
    description = 'Select sequence segments in the middle column for beta and visual arrestins. You can expand every structural element and select individual' \
        + ' residues by clicking on the down arrows next to each helix, sheet or loop.\n\n You can select the full sequence or show all structured regions at the same time.\n\nSelected segments will appear in the' \
        + ' right column, where you can edit the list.\n\nOnce you have selected all your segments, click the green' \
        + ' button.'

    template_name = 'common/segmentselection.html'

    selection_boxes = OrderedDict([
        ('reference', False),
        ('targets', True),
        ('segments', True),
    ])
    buttons = {
        'continue': {
            'label': 'Show alignment',
            'url': '/alignment/render',
            'color': 'success',
        },
    }

    position_type = 'arrestin'

    ## Add some Arrestin specific positions
    rsets = ResiduePositionSet.objects.filter(name__in=['Arrestin interface']).prefetch_related('residue_position')

    ## ProteinSegment for different proteins
    ss = ProteinSegment.objects.filter(partial=False, proteinfamily='Arrestin').prefetch_related('generic_numbers')
    ss_cats = ss.values_list('category').order_by('category').distinct('category')


class BlastSearchInput(AbsMiscSelection):
    step = 1
    number_of_steps = 1
    docs = 'sequences.html#similarity-search-blast'
    title = 'BLAST search'
    description = 'Enter a sequence into the text box and press the green button.'
    buttons = {
        'continue': {
            'label': 'BLAST',
            'onclick': 'document.getElementById(\'form\').submit()',
            'color': 'success',
        },
    }
    selection_boxes = {}
    blast_input = True


class BlastSearchResults(TemplateView):
    """
    An interface for blast similarity search of the input sequence.
    """
    template_name="blast/blast_search_results.html"

    def post(self, request, *args, **kwargs):

        if 'human' in request.POST.keys():
            blast = BlastSearch(blastdb=os.sep.join([settings.STATICFILES_DIRS[0], 'blast', 'protwis_human_blastdb']), top_results=50)
            blast_out = blast.run(request.POST['input_seq'])
        else:
            blast = BlastSearch(top_results=50)
            blast_out = blast.run(request.POST['input_seq'])

        context = {}
        context['results'] = [(Protein.objects.get(pk=x[0]), x[1]) for x in blast_out]
        context["input"] = request.POST['input_seq']

        return render(request, self.template_name, context)


def render_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)
    if simple_selection == False or not simple_selection.targets:
        return redirect("/alignment/targetselection")

    # create an alignment object
    a = Alignment()

    # load data from selection into the alignment
    # only show wildtype protein entries if selection type is family
    if len([t for t in simple_selection.targets if t.type=='family'])>0:
        a.load_proteins_from_selection(simple_selection, only_wildtype=True)
    else:
        a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    key = "ALIGNMENT_" + a.get_hash()
    return_html = cache_alignment.get(key)

    if return_html==None:
        # build the alignment data matrix
        check = a.build_alignment()
        if check == 'Too large':
            return render(request, 'alignment/error.html', {'proteins': len(a.proteins), 'residues':a.number_of_residues_total})

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()
        num_of_sequences = len(a.proteins)
        num_residue_columns = len(a.positions) + len(a.segments)

        return_html = render(request, 'alignment/alignment.html', {'a': a, 'num_of_sequences': num_of_sequences,
            'num_residue_columns': num_residue_columns})

    cache_alignment.set(key, return_html, 60*60*24*7) #set alignment cache one week

    return return_html

def render_family_alignment(request, slug):
    # create an alignment object
    a = Alignment()

    # fetch proteins and segments
    proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt')

    if len(proteins)>50 and len(slug.split("_"))<4:
        # If alignment is going to be too big, only pick human.
        proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt', species__latin_name='Homo sapiens')

    if slug.startswith('100'):

        gsegments = definitions.G_PROTEIN_SEGMENTS

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(gsegments['Full'])])
        segments = ProteinSegment.objects.filter(slug__in=gsegments['Full'], partial=False).order_by(preserved)

    elif slug.startswith('200'):
        arrsegments = definitions.ARRESTIN_SEGMENTS

        preserved = Case(*[When(slug=pk, then=pos) for pos, pk in enumerate(arrsegments['Full'])])
        segments = ProteinSegment.objects.filter(slug__in=arrsegments['Full'], partial=False).order_by(preserved)

    else:
        segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR')
        if len(proteins)>50:
            # if a lot of proteins, exclude some segments
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(slug__in=['N-term','C-term'])
        if len(proteins)>200:
            # if many more proteins exluclude more segments
            segments = ProteinSegment.objects.filter(partial=False, proteinfamily='GPCR').exclude(slug__in=['N-term','C-term']).exclude(category='loop')

    protein_ids = []
    for p in proteins:
        protein_ids.append(p.pk)
    protein_list = ','.join(str(x) for x in sorted(protein_ids))

    #create unique proteins_id
    segments_ids = []
    for s in segments:
        segments_ids.append(s.slug)
    segments_list = ','.join(str(x) for x in sorted(segments_ids))

    # Store proteins and segments as selection to enable Fasta/Excel/CSV downloads
    selection = Selection()
    for prot in proteins:
        selection.add('targets', 'protein', SelectionItem('protein', prot))
    for segment in segments:
        selection.add('segments', 'protein_segment', SelectionItem('protein_segment', segment))
    request.session['selection'] = selection.exporter()

    s = str(protein_list+"_"+segments_list)
    key = "ALIGNMENT_"+hashlib.md5(s.encode('utf-8')).hexdigest()
    return_html = cache_alignment.get(key)
    if return_html==None:
        # load data into the alignment
        a.load_proteins(proteins)
        a.load_segments(segments)

        # build the alignment data matrix
        a.build_alignment()

        # calculate consensus sequence + amino acid and feature frequency
        a.calculate_statistics()

        num_of_sequences = len(a.proteins)
        num_residue_columns = len(a.positions) + len(a.segments)

        return_html = render(request, 'alignment/alignment.html', {'a': a, 'num_of_sequences': num_of_sequences,
        'num_residue_columns': num_residue_columns})

    #update it if used
    cache_alignment.set(key,return_html, 60*60*24*7) #set alignment cache one week

    return return_html

def render_fasta_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    response = render(request, 'alignment/alignment_fasta.html', context={'a': a}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.fasta"
    return response

def render_fasta_family_alignment(request, slug):
    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # fetch proteins and segments
    proteins = Protein.objects.filter(family__slug__startswith=slug, sequence_type__slug='wt')
    segments = ProteinSegment.objects.filter(partial=False)

    # load data into the alignment
    a.load_proteins(proteins)
    a.load_segments(segments)

    # build the alignment data matrix
    a.build_alignment()

    response = render(request, 'alignment/alignment_fasta.html', context={'a': a}, content_type='text/fasta')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.fasta"
    return response

def render_csv_alignment(request):
    # get the user selection from session
    simple_selection = request.session.get('selection', False)

    # create an alignment object
    a = Alignment()
    a.show_padding = False

    # load data from selection into the alignment
    a.load_proteins_from_selection(simple_selection)
    a.load_segments_from_selection(simple_selection)

    # build the alignment data matrix
    a.build_alignment()

    # calculate consensus sequence + amino acid and feature frequency
    a.calculate_statistics()

    response = render(request, 'alignment/alignment_csv.html', context={'a': a}, content_type='text/csv')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.csv"
    return response

# Excel download based on seq. signature tool
def render_alignment_excel(request):

    # Grab all targets
    targets = request.session.get('selection', False)

    # create placeholder seq signature
    signature = SequenceSignature()
    signature.setup_alignments_from_selection(targets, targets)

    # calculate the signture
    signature.calculate_signature()
    signature.calculate_zscales_signature()

    outstream = BytesIO()
    wb = xlsxwriter.Workbook(outstream, {'in_memory': True})

    # Sequence alignment of targets
    signature.prepare_excel_worksheet(
        wb,
        'Alignment',
        'positive',
        'alignment'
    )

    # Residue properties stats
    signature.prepare_excel_worksheet(
        wb,
        'Property_conservation',
        'positive',
        'features'
    )
    # Z-scales
    signature.zscales_excel(
        wb,
        "Z-scales",
        'positive'
    )
    wb.close()
    outstream.seek(0)
    response = HttpResponse(
        outstream.read(),
        content_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        )
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_alignment.xlsx"

    return response

def retrieve_class_similarity_matrix(output_type='html',classless=False,human_only=False):
    cross_class_similarities = OrderedDict()
    class_representative_species_list = ClassRepresentativeSpecies.objects.all().prefetch_related()
    if human_only:
        class_representative_species_list = class_representative_species_list.filter(species__common_name__iexact='human')
        human_classes_set = set()
        for class_representative_species in class_representative_species_list:
            human_classes_set.add(class_representative_species.protein_family)
    if output_type == 'html':
        gpcr_family_name_2_class_representative_species = {}
        for class_representative_species in class_representative_species_list:
            gpcr_family_name = class_representative_species.protein_family.name
            species_name = class_representative_species.species.common_name
            if species_name.lower() != 'human':
                if strain_re.search(species_name):
                    species_name = '<i>'+escape(class_representative_species.species.latin_name)+'</i> <br>'+escape(species_name)
                    gpcr_family_name_2_class_representative_species[gpcr_family_name] = species_name
                else:
                    species_name = class_representative_species.species.latin_name
                    gpcr_family_name_2_class_representative_species[gpcr_family_name] = '<i>'+escape(species_name)+'</i>'
            else:
                gpcr_family_name_2_class_representative_species[gpcr_family_name] = escape(species_name)

    ties_list = ClassSimilarityTie.objects.all().select_related('protein1','protein2').values('class_similarity_id',
                                                                    'protein1__entry_name','protein2__entry_name',
                                                                    'protein1__name','protein2__name','type')

    class_similarity_id_2_ties = {}
    for tie in ties_list:
        class_similarity_id = tie['class_similarity_id']
        if class_similarity_id not in class_similarity_id_2_ties:
            class_similarity_id_2_ties[class_similarity_id] = []
        class_similarity_id_2_ties[class_similarity_id].append(tie)

    class_sim = ClassSimilarity.objects.all().select_related('protein_family1','protein_family2')
    if human_only:
        class_sim = class_sim.filter(protein_family1__in=human_classes_set,protein_family2__in=human_classes_set)


    if not classless:
        for slug in list(CLASSLESS_PARENT_GPCR_SLUGS):
            class_sim = class_sim.exclude(protein_family1__slug__startswith=slug)
            class_sim = class_sim.exclude(protein_family2__slug__startswith=slug)


    for class_pair in class_sim.values('id','protein_family1__slug','protein_family1__name',
                                                                   'protein_family2__slug','protein_family2__name',
                                                                    'similar_protein1__entry_name','similar_protein2__entry_name',
                                                                    'similar_protein1__name','similar_protein2__name',
                                                                    'ident_protein1__entry_name','ident_protein2__entry_name',
                                                                    'ident_protein1__name','ident_protein2__name','identity','similarity')\
                                                                    .order_by('protein_family1__slug','protein_family2__slug'):
        gpcr_class1_name = class_pair['protein_family1__name']
        gpcr_class2_name = class_pair['protein_family2__name']
        if gpcr_class1_name not in cross_class_similarities:
            cross_class_similarities[gpcr_class1_name] = OrderedDict()
        if gpcr_class2_name not in cross_class_similarities[gpcr_class1_name]:
            cross_class_similarities[gpcr_class1_name][gpcr_class2_name] = {}
        cross_class_sim = cross_class_similarities[gpcr_class1_name][gpcr_class2_name]
        cross_class_sim['identity'] = class_pair['identity']
        cross_class_sim['similarity'] = class_pair['similarity']
        if output_type == 'html':
            name_type = 'name'
            ident_protein1_name = Protein(name=class_pair['ident_protein1__'+name_type]).short().replace('<i>','').replace('</i>','')
            ident_protein2_name = Protein(name=class_pair['ident_protein2__'+name_type]).short().replace('<i>','').replace('</i>','')
            similar_protein1_name = Protein(name=class_pair['similar_protein1__'+name_type]).short().replace('<i>','').replace('</i>','')
            similar_protein2_name = Protein(name=class_pair['similar_protein2__'+name_type]).short().replace('<i>','').replace('</i>','')
            name_type = 'entry_name'
            ident_protein1_entry_name = class_pair['ident_protein1__'+name_type]
            ident_protein2_entry_name = class_pair['ident_protein2__'+name_type]
            similar_protein1_entry_name = class_pair['similar_protein1__'+name_type]
            similar_protein2_entry_name = class_pair['similar_protein2__'+name_type]
            cross_class_sim['identity_gpcr_pair_entry_name'] = (ident_protein1_entry_name,ident_protein2_entry_name)
            cross_class_sim['identity_gpcr_pair_w_entry_name'] = []
            cross_class_sim['similarity_gpcr_pair'] = (similar_protein1_name,similar_protein2_name)
            cross_class_sim['similarity_gpcr_pair_entry_name'] = (similar_protein1_entry_name,similar_protein2_entry_name)
            cross_class_sim['similarity_gpcr_pair_w_entry_name'] = []
        else:
            name_type = 'entry_name'
            ident_protein1_name = class_pair['ident_protein1__'+name_type]
            ident_protein2_name = class_pair['ident_protein2__'+name_type]
            similar_protein1_name = class_pair['similar_protein1__'+name_type]
            similar_protein2_name = class_pair['similar_protein2__'+name_type]
        cross_class_sim['identity_gpcr_pair'] = (ident_protein1_name,ident_protein2_name)
        cross_class_sim['identity_gpcr_pair_w'] = []
        cross_class_sim['similarity_gpcr_pair'] = (similar_protein1_name,similar_protein2_name)
        cross_class_sim['similarity_gpcr_pair_w'] = []


        if class_pair['id'] in class_similarity_id_2_ties:
            for tie in class_similarity_id_2_ties[class_pair['id']]:
                if tie['type'] == ClassSimilarityType.IDENTITY:
                    type = 'identity'
                elif tie['type'] == ClassSimilarityType.SIMILARITY:
                    type = 'similarity'
                
                if output_type == 'html':
                    name_type = 'name'
                    protein1 = Protein(name=tie['protein1__'+name_type]).short().replace('<i>','').replace('</i>','')
                    protein2 = Protein(name=tie['protein2__'+name_type]).short().replace('<i>','').replace('</i>','')
                    name_type = 'entry_name'
                    protein1_entry_name= tie['protein1__'+name_type]
                    protein2_entry_name = tie['protein2__'+name_type]
                    cross_class_sim[type+'_gpcr_pair_w_entry_name'].append((protein1_entry_name,protein2_entry_name))
                else:
                    name_type = 'entry_name'
                    protein1 = tie['protein1__'+name_type]
                    protein2 = tie['protein2__'+name_type]

                cross_class_sim[type+'_gpcr_pair_w'].append((protein1,protein2))

    selected_parent_gpcr_families_names1 = OrderedDict()
    selected_parent_gpcr_families_names2 = OrderedDict()
    for key in cross_class_similarities:
        selected_parent_gpcr_families_names1[key] = None
        for key2 in cross_class_similarities[key]:
            selected_parent_gpcr_families_names2[key2] = None
    for key in selected_parent_gpcr_families_names2:
        selected_parent_gpcr_families_names1[key] = None

    selected_parent_gpcr_families_names = selected_parent_gpcr_families_names1.keys()
    
    cross_class_similarity_matrix = OrderedDict()
    i = 0
    for key in selected_parent_gpcr_families_names:
        row = []
        j = 0
        for key2 in selected_parent_gpcr_families_names:
            if key == key2:
                row.append(['-','-',''])
                continue
            if i > j:
                type = 'similarity'
            else:
                type = 'identity'
            if key in cross_class_similarities:
                if key2 in cross_class_similarities[key]:
                    sim = cross_class_similarities[key][key2]
                else:
                    sim = cross_class_similarities[key2][key]
            else:
                sim = cross_class_similarities[key2][key]
            pairs = [sim[type+'_gpcr_pair']]
            pairsw = sim[type+'_gpcr_pair_w']
            if len(pairsw) > 0:
                pairs += pairsw
            if output_type == 'html':
                pairs_entry_name = [sim[type+'_gpcr_pair_entry_name']]
                pairsw_entry_name = sim[type+'_gpcr_pair_w_entry_name']
                if len(pairsw_entry_name) > 0:
                    pairs_entry_name += pairsw_entry_name
                row.append([str(sim[type]),str(sim[type]//10),pairs,pairs_entry_name])
                species_name = gpcr_family_name_2_class_representative_species[key] 
            else:
                row.append([str(sim[type]),str(sim[type]//10),pairs])
            j += 1
        
        name = key.replace('<i>','').replace('</i>','')
        if output_type == 'html':
            name = class_fullname_re.sub(r'\1<br>\3', name)
            name = class_fungal_re.sub(r'\1<br>\3<br>\5', name)
            cross_class_similarity_matrix[key] = {'name':name,'values':row, 'species':species_name}
        else:
            cross_class_similarity_matrix[key] = {'name':name,'values':row}
        i += 1 
    return (cross_class_similarity_matrix,selected_parent_gpcr_families_names)

render_class_similarity_csv_matrix_urls = ['human_only_without_classless', 'human_only_with_classless','all_without_classless','all_with_classless']

def render_class_similarity_matrix(request):
    r_human_only_without_classless = retrieve_class_similarity_matrix(classless=False,human_only=True)
    r_human_only_with_classless = retrieve_class_similarity_matrix(classless=True,human_only=True)
    r_all_without_classless = retrieve_class_similarity_matrix(classless=False,human_only=False)
    r_all_with_classless = retrieve_class_similarity_matrix(classless=True,human_only=False)

    list_tmp = [r_human_only_without_classless,r_human_only_with_classless,r_all_without_classless,r_all_with_classless]

    h_list = ['Human only', 'Human only (including Classless)','All (including non-human)','All (including non-human & Classless)']
    csv_list = [{'classless': '0','human_only':'1'},{'classless': '1','human_only':'1'},{'classless': '0','human_only':'0'},{'classless': '1','human_only':'0'}]

    return render(request, 'class_similarity/matrix.html', {'csv_list':csv_list,'h_list': h_list, 'p_list': [r[1] for r in list_tmp],'m_list': [r[0] for r in list_tmp]})



def render_class_similarity_csv_matrix(request):
    classless = True
    human_only = False
    try:
        classless_g = int(request.GET.get('classless', 1))
        if classless_g == 0:
            classless = False
    except ValueError as e:
        pass
    try:
        human_only_g = int(request.GET.get('human_only', 0))
        if abs(human_only_g) > 0:
            human_only = True
    except ValueError as e:
        pass

    r = retrieve_class_similarity_matrix(output_type='csv',classless=classless,human_only=human_only)

    response = render(request, 'class_similarity/matrix_csv.html', {'p': r[1],'m': r[0]}, content_type='text/csv')
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_similaritymatrix.csv"
    return response

def render_class_similarity_xlsx_matrix(request):
    classless = True
    human_only = False
    try:
        classless_g = int(request.GET.get('classless', 1))
        if classless_g == 0:
            classless = False
    except ValueError as e:
        pass
    try:
        human_only_g = int(request.GET.get('human_only', 0))
        if abs(human_only_g) > 0:
            human_only = True
    except ValueError as e:
        pass

    r = retrieve_class_similarity_matrix(output_type='csv',classless=classless,human_only=human_only)
    m, p = r

    xlsx_output = BytesIO()
    
    workbook = xlsxwriter.Workbook(xlsx_output, {'in_memory': True})

    axis = True

    v_axis_row = 1    #starting row index of vertical axis
    v_axis_col = 0    #column index of identity axis
    h_axis_row = 0    #row index of horitzontal axis
    h_axis_col = 1    #starting column index of indentity axis
    row_header_row = v_axis_row + 1    #starting row index of the row header
    row_header_col = v_axis_col + 1    #column index of the row header
    col_header_row = h_axis_row + 1    #starting column index of the column header
    col_header_col = h_axis_col + 1    #row index of the column header
    if not axis:
        row_header_row = 1    #starting row index of the row header
        row_header_col = 0    #column index of the row header
        col_header_row = 0    #starting column index of the column header
        col_header_col = 1    #row index of the column header
    
    
    first_data_row = row_header_row    #first data cell row index
    first_data_col = col_header_col    #first data cell column index

    n_header = len(m.keys())

    last_data_row = n_header - 1 + first_data_row    #last data cell row index
    last_data_col = n_header - 1 + first_data_col    #last data cell column index

    #format

    table_headers = [v['name'] for p, v in m.items()]
    worksheet1 = workbook.add_worksheet('Values')
    worksheet2 = workbook.add_worksheet('Receptors')

    table_headers_format = workbook.add_format()            #row and column headers format
    table_h_axis_format = workbook.add_format()             #horitzontal axis format
    table_v_axis_format = workbook.add_format()             #vertical axis format
    table_headers_numbers_format = workbook.add_format()    #row and column headers format for cells containing numeric data
    numbers_format = workbook.add_format()                  #data cell format for cells containing numeric data                                        
    cell_format = workbook.add_format()                     #data cell format for cells containing non-numeric data 

    table_headers_format.set_bold(True)
    table_headers_numbers_format.set_bold(True)
    table_h_axis_format.set_bold(True)
    table_v_axis_format.set_bold(True)

    table_headers_format.set_align('vcenter')
    table_headers_numbers_format.set_align('vcenter')
    table_h_axis_format.set_align('vcenter')
    table_v_axis_format.set_align('vcenter')
    numbers_format.set_align('vcenter')
    cell_format.set_align('vcenter')

    table_headers_numbers_format.set_align('center')
    table_h_axis_format.set_align('center')
    table_v_axis_format.set_align('center')
    numbers_format.set_align('center')

    numbers_format.set_border(1)
    numbers_format.set_border_color('#DDDDDD')
    table_v_axis_format.set_rotation(90)
    table_headers_format.set_text_wrap()
    table_headers_numbers_format.set_text_wrap()

    #axis titles

    if axis:
            worksheet1.merge_range(h_axis_row , h_axis_col, h_axis_row, last_data_col, 'Identity (%)',table_h_axis_format)
            worksheet2.merge_range(h_axis_row , h_axis_col, h_axis_row, last_data_col, 'Identity',table_h_axis_format)
    if axis:
            worksheet1.merge_range(v_axis_row , v_axis_col, last_data_row, v_axis_col, 'Similarity (%)',table_v_axis_format)
            worksheet2.merge_range(v_axis_row , v_axis_col, last_data_row, v_axis_col, 'Similarity',table_v_axis_format)

    #column header
    worksheet1.write_row(col_header_row , col_header_col, table_headers, table_headers_numbers_format)
    worksheet2.write_row(col_header_row, col_header_col, table_headers, table_headers_format)

    row = first_data_row
    for p, v in m.items():
        worksheet1.write(row,row_header_col,v['name'],table_headers_numbers_format)
        col = first_data_col 
        for v_col in v['values']:
            cell_value = v_col[0]
            if cell_value  != '-':
                cell_value = int(cell_value )
            worksheet1.write(row,col,cell_value,numbers_format)
            col += 1
        row += 1


    row = first_data_row
    for p, v in m.items():
        worksheet2.write(row,row_header_col,v['name'],table_headers_format)
        col = first_data_col
        for v_col in v['values']:
            cell_value = ''
            for counter, pair in enumerate(v_col[2]):
                if row > col:
                    cell_value += pair[1]+' vs '+pair[0]
                elif row < col:
                    cell_value += pair[0]+' vs '+pair[1]
                if counter != len(v_col[2]) - 1:
                    cell_value += ':'
            if row == col:
                cell_value = '-'
            worksheet2.write(row,col,cell_value,cell_format)
            col += 1
        row += 1
    try:
        # Requires xlsxwriter > 3.0.8

        worksheet1.autofit()
        worksheet2.autofit()
    except Exception as e:
        if axis:
            first_col = v_axis_col
        else:
            first_col = row_header_col
        worksheet2.set_column(first_col, last_data_col, 12)
        worksheet2.set_column(first_data_col, last_data_col, 81)
        pass
    worksheet1.set_column(0, last_data_col, 12)
    for row in range(0,last_data_col + 1):
        worksheet1.set_row(row, 50)
    worksheet1.freeze_panes(2, 2)
    worksheet2.freeze_panes(2, 2)
    multi_range_list = []
    # compute cell multi-range for color scale

    # identity
    multi_range_list = []
    for i in range(1,last_data_row - first_data_row + 1):
        first_row = first_data_row
        last_row = first_data_row + i - 1 # do not color scale the diagonal
        first_col = last_col = first_data_col + i
        multi_range_list.append(xl_range_abs(first_row,first_col,last_row,last_col))
    worksheet1.conditional_format(multi_range_list[0],
                                {'type': '2_color_scale','min_color': '#FFFFFF',
                                'max_color': '#999999','multi_range': ' '.join(multi_range_list)})

    # similarity
    multi_range_list = []
    for i in range(0,last_data_row - first_data_row + 1):
        first_row = first_data_row + 1 + i # do not color scale the diagonal
        last_row = last_data_row
        first_col = last_col = first_data_col + i
        multi_range_list.append(xl_range_abs(first_row,first_col,last_row,last_col))
    worksheet1.conditional_format(multi_range_list[0],
                                  {'type': '2_color_scale','min_color': '#FFFFFF',
                                    'max_color': '#999999','multi_range': ' '.join(multi_range_list)})
    workbook.close()
    xlsx_file = File(xlsx_output)
    xlsx_file_size = xlsx_file.size
    
    response = HttpResponse(xlsx_file ,'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    response['Content-Length'] = xlsx_file_size
    response['Content-Disposition'] = "attachment; filename=" + site_title(request)["site_title"] + "_similaritymatrix.xlsx"
    return response

