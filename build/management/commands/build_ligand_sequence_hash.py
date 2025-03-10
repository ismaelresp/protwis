from django.core.management.base import BaseCommand, CommandError
from django.db import connection
from django.db import IntegrityError
from django.db.models import Count, Max, F, Subquery, OuterRef

from psycopg2.errors import UniqueViolation

from ligand.models import Ligand

import os
import django.apps
import logging
import hashlib
import base64
from collections import OrderedDict

max_buffer_size = 1000

class Command(BaseCommand):
    help = 'Build ligand sequence hash. '

    def add_arguments(self, parser):
        super(Command, self).add_arguments(parser=parser)
        parser.add_argument('--verbose', default=False, action='store_true', help='Print progress in stdout.')

    logger = logging.getLogger(__name__)


    def handle(self, *args, **options):    
        error = None
        if options['verbose']: print('Building ligand sequence hashes...')

        query_fields = ['id','sequence']
        update_fields = ['sequence_hash','sequence_hash_col','sequence_dup','sequence_hash_and_col_main']
        q = Ligand.objects.exclude(sequence=None).order_by('id').values(*query_fields)
        i = 0
        while True:
            q_results = list(q[i:max_buffer_size+i])
            print('hola:',len(q_results))
            if not q_results:
                break
            objs = []
            for q_result in q_results:    
                hash = base64.b32encode(hashlib.md5(q_result['sequence'].encode()).digest()).decode().strip('=')
                q_result['sequence_hash'] = hash
                q_result['sequence_hash_col'] = ''
                q_result['sequence_dup'] = 0
                q_result['sequence_hash_and_col_main'] = True
                objs.append(Ligand(**q_result))
            # Try to save in DB the hashes
            try:
                update_q = Ligand.objects.bulk_update(objs,update_fields)
            except (IntegrityError,UniqueViolation) as e:


                # Separate in two different dict() Ligand objects with duplicated hashes/collisions from the other hashes 
                objs_dict = {j:obj for j,obj in enumerate(objs)}
                del objs
                sequence_hashes_dict = {}
                dup_sequence_hashes_set = set()
                for j,q_result in enumerate(q_results):
                    hash = q_result['sequence_hash']
                    if hash in sequence_hashes_dict:
                        dup_sequence_hashes_set.add(hash)
                        if hash not in sequence_hashes_dict:
                            sequence_hashes_dict[hash] = []
                        sequence_hashes_dict[hash].append(objs_dict[j])
                    else:
                        sequence_hashes_dict[hash] = [objs_dict[j]]
                del objs_dict

                print('hola2:','2SLRRWX6PFJMBE46AQVIJHQEVI' in dup_sequence_hashes_set)

                list_of_unique_hashes = list(sequence_hashes_dict.keys())

                q_max_sequence_dup = Ligand.objects.filter(sequence_hash__in=list_of_unique_hashes)
                q_max_sequence_dup = q_max_sequence_dup.values('sequence_hash','sequence_hash_col')
                q_max_sequence_dup = q_max_sequence_dup.annotate(max_sequence_dup=Max('sequence_dup'))
                # q_max_sequence_dup = q_max_sequence_dup.filter(sequence_dup=F('max_sequence_dup'))
                q_max_sequence_dup_sub = q_max_sequence_dup.filter(sequence_hash=OuterRef('sequence_hash'))
                q_max_sequence_dup_sub = q_max_sequence_dup_sub.filter(sequence_hash_col=OuterRef('sequence_hash_col'))
                q_max_sequence_dup_sub = q_max_sequence_dup_sub.values('max_sequence_dup')
                # q_num_col_or_dup = q_max_sequence_dup.annotate(num_col_or_dup=Count('id'))
                

                q_num_col_or_dup = Ligand.objects.filter(sequence_dup=Subquery(q_max_sequence_dup_sub))
                q_num_col_or_dup = q_num_col_or_dup.values('sequence_hash').annotate(num_col_or_dup=Count('id'))

                # q_num_col_or_dup = Ligand.objects.filter(sequence_hash__in=list_of_unique_hashes,
                #                                         sequence_hash=q_max_sequence_dup['sequence_hash'],
                #                                         sequence_hash_col=q_max_sequence_dup['sequence_hash_col'],
                #                                         sequence_dup=q_max_sequence_dup['max_sequence_dup'],
                #                                         )
                # q_num_col_or_dup = q_num_col_or_dup.values('sequence_hash')
                # q_num_col_or_dup = q_num_col_or_dup.values('sequence_hash').annotate(num_col_or_dup=Count('id'))
                
                # Make a list of batches of sequence hashes than will return a query with a number of records < max_buffer_size
                # The following code does not preserve the order of the sequence hashes
                col_or_dup_batch_list = []
                overrun_buffer_hashes_list = []
                current_batch_size = 0
                batch = []
                for row in q_num_col_or_dup:
                    num_col_or_dup = row['num_col_or_dup']
                    if num_col_or_dup > max_buffer_size:
                        sequence_hash = row['sequence_hash']
                        # Remove this warning if the code already takes care of this
                        self.logger.warning("build_ligand_sequence_hash: Number of ligand duplicates + number of " + \
                                            "hash collisions larger than %d for %s." % (max_buffer_size, sequence_hash) + \
                                            "This might cause RAM memory overrun during building with the current " + \
                                            "implementation." )
                        overrun_buffer_hashes_list.append(sequence_hash)
                        continue
                    
                    current_batch_size += num_col_or_dup
                    if current_batch_size > max_buffer_size:
                        col_or_dup_batch_list.append(batch)
                        batch = []
                        current_batch_size = num_col_or_dup
                    batch.append(row['sequence_hash'])
                if len(batch) > 0:
                    col_or_dup_batch_list.append(batch)

                # Today, we still don't take care of buffer overruning hashes
                col_or_dup_batch_list += overrun_buffer_hashes_list
                del overrun_buffer_hashes_list

                fixed_dup_or_col_hashes = set()
                for batch in col_or_dup_batch_list:
                    
                    q_max_sequence_dup = Ligand.objects.filter(sequence_hash__in=batch)
                    q_max_sequence_dup = q_max_sequence_dup.values('sequence_hash','sequence_hash_col')
                    q_max_sequence_dup = q_max_sequence_dup.annotate(max_sequence_dup=Max('sequence_dup'))
                    # q_max_sequence_dup = q_max_sequence_dup.filter(sequence_dup=F('max_sequence_dup'))
                    q_max_sequence_dup_sub = q_max_sequence_dup.filter(sequence_hash=OuterRef('sequence_hash'))
                    q_max_sequence_dup_sub = q_max_sequence_dup_sub.filter(sequence_hash_col=OuterRef('sequence_hash_col'))
                    q_max_sequence_dup_sub = q_max_sequence_dup_sub.values('max_sequence_dup')
                                                        

                    q_col_or_dup = Ligand.objects.filter(sequence_dup=Subquery(q_max_sequence_dup_sub))
                    q_col_or_dup = q_col_or_dup.values(*query_fields+update_fields)

                    q_col_or_dup_hash_seq_dict = {}
                    for r in q_col_or_dup:
                        hash = r['sequence_hash']
                        if hash not in q_col_or_dup_hash_seq_dict:
                            q_col_or_dup_hash_seq_dict[hash] = {}
                        q_col_or_dup_hash_seq_dict[hash][r['sequence']] = r
                    del q_col_or_dup

                    for hash, seq_dict in q_col_or_dup_hash_seq_dict.items():
                        # hash_cols is an hexadecimal number written from right to left
                        hash_cols_decimal_max = max([int(r['sequence_hash_col'][::-1],16) for r in seq_dict.values()\
                                                      if r['sequence_hash_col'] != '']+[0])
                        
                            
                        db_col_dup_count_dict = {}
                        new_col_dup_count_dict = {j:0 for j in range(0,hash_cols_decimal_max)}
                        new_col_count = 0
                        new_col_sequences_dict = {}
                        for obj in sequence_hashes_dict[hash]:
                            if '2SLRRWX6PFJMBE46AQVIJHQEVI' == hash:
                                print('hola3:','2SLRRWX6PFJMBE46AQVIJHQEVI' == hash)
                            # check if it is a duplicate
                            if  obj.sequence in seq_dict:
                                if '2SLRRWX6PFJMBE46AQVIJHQEVI' == hash:
                                    print('hola4:','2SLRRWX6PFJMBE46AQVIJHQEVI' == hash)
                                r = seq_dict[obj.sequence]
                                hash_col = r['sequence_hash_col']
                                obj.sequence_hash_col = hash_col
                                if hash_col not in db_col_dup_count_dict:
                                    db_col_dup_count_dict[hash_col] = r['sequence_dup']
                                db_col_dup_count_dict[hash_col] += 1
                                obj.sequence_dup = db_col_dup_count_dict[hash_col]
                                obj.sequence_hash_and_col_main = False
                                


                            elif obj.sequence in new_col_sequences_dict:
                                if '2SLRRWX6PFJMBE46AQVIJHQEVI' == hash:
                                    print('hola5:','2SLRRWX6PFJMBE46AQVIJHQEVI' == hash)
                                hash_col = new_col_sequences_dict[obj.sequence]
                                obj.sequence_hash_col = hash_col
                                new_col_dup_count_dict[new_col_count] += 1
                                obj.sequence_dup = new_col_dup_count_dict[new_col_count]
                                obj.sequence_hash_and_col_main = False
                            else:
                                if '2SLRRWX6PFJMBE46AQVIJHQEVI' == hash:
                                    print('hola6:','2SLRRWX6PFJMBE46AQVIJHQEVI' == hash)
                                new_col_count += 1
                                hash_col = hex(hash_cols_decimal_max + new_col_count)[2:][::-1]
                                obj.sequence_hash_col = hash_col
                                obj.sequence_hash_and_col_main = True
                                obj.sequence_dup = 0
                                new_col_sequences_dict[obj.sequence] = hash_col
                        if hash_cols_decimal_max == 0:
                            if '2SLRRWX6PFJMBE46AQVIJHQEVI' == hash:
                                print('hola7:','2SLRRWX6PFJMBE46AQVIJHQEVI' == hash)
                            sequences = list(seq_dict.keys())
                            if len(sequences) > 1:
                                self.logger.warning("build_ligand_sequence_hash: sequence hash %s " % (sequence_hash) + \
                                            "had collisions and sequence_hash_col field is empty.")
                            r = seq_dict[sequences[0]]
                            obj = Ligand(**r)
                            obj.sequence_hash_col = '0'
                            sequence_hashes_dict[hash].append(obj)
                        fixed_dup_or_col_hashes.add(hash)

                dup_or_col_hashes_not_in_db = dup_sequence_hashes_set - fixed_dup_or_col_hashes  

                # Assign duplicated hashes/collisions IDs
                for hash in list(dup_or_col_hashes_not_in_db):
                    dup_col_objs = sequence_hashes_dict[hash]
                    new_col_dup_count_dict = {0:0}
                    new_col_count = 0
                    new_col_sequences_dict = {}
                    for obj in dup_col_objs:
                        if obj.sequence in new_col_sequences_dict:
                            hash_col = new_col_sequences_dict[obj.sequence]
                            obj.sequence_hash_col = hash_col
                            new_col_dup_count_dict[new_col_count] += 1
                            obj.sequence_dup = new_col_dup_count_dict[new_col_count]
                            obj.sequence_hash_and_col_main = False
                        else:
                            hash_col = hex(new_col_count)[2:][::-1]
                            obj.sequence_hash_col = hash_col
                            obj.sequence_hash_and_col_main = True
                            obj.sequence_dup = 0
                            new_col_sequences_dict[obj.sequence] = hash_col
                            new_col_count += 1
                            new_col_dup_count_dict[new_col_count] = 0
                    if len(new_col_sequences_dict.keys()) < 2:
                        for obj in dup_col_objs:
                            obj.sequence_hash_col = ''

                # Convert into a list the dict() with a list of Ligand objects with duplicated hashes/collisions
                objs = []
                for v in sequence_hashes_dict.values():
                    objs += v
                del sequence_hashes_dict
                update_q = Ligand.objects.bulk_update(objs,update_fields)
                pass


            i += max_buffer_size
            
            
        self.logger.info('Ligand sequence hashes built.')
            



 