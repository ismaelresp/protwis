from django.db import models
from django.db import connection
from enum import IntEnum


# Create your models here.

class AlignmentConsensus(models.Model):
    slug = models.SlugField(max_length=100, unique=True)
    alignment = models.BinaryField()
    gn_consensus = models.BinaryField(blank=True) # Store conservation calculation for each GN


class CustomClassSimilarityManager(models.Manager):
    def truncate_table(self):
        cursor = connection.cursor()
        table_name = self.model._meta.db_table
        sql = 'TRUNCATE TABLE "{0}" CASCADE'.format(table_name)
        cursor.execute(sql)

class ClassSimilarity(models.Model):
    protein_family1 = models.ForeignKey('protein.ProteinFamily',null=False, on_delete=models.CASCADE, related_name='class_similarity_protein_family1')
    protein_family2 = models.ForeignKey('protein.ProteinFamily',null=False, on_delete=models.CASCADE, related_name='class_similarity_protein_family2')
    similar_protein1 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_similar_protein1')
    similar_protein2 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_similar_protein2')
    ident_protein1 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_ident_protein1')
    ident_protein2 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_ident_protein2')
    similarity = models.IntegerField(null=False)
    identity = models.IntegerField(null=False)

    objects = models.Manager()  # The default manager.
    custom_objects = CustomClassSimilarityManager()  # The custom manager.
    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['protein_family1', 'protein_family2'], name='unique_class_similarity_protein_family'),
            models.UniqueConstraint(fields=['similar_protein1', 'similar_protein2'], name='unique_class_similarity_similar_proteins'),
            models.UniqueConstraint(fields=['ident_protein1', 'ident_protein2'], name='unique_class_similarity_ident_proteins'),
        ]
    def __str__(self):
        return str(self.protein_family1)+" - "+str(self.protein_family2)+": "+str(self.similarity)

class CustomClassRepresentativeSpeciesManager(models.Manager):
    def truncate_table(self):
        cursor = connection.cursor()
        table_name = self.model._meta.db_table
        sql = 'TRUNCATE TABLE "{0}"'.format(table_name)
        cursor.execute(sql)

class ClassRepresentativeSpecies(models.Model):
    protein_family = models.OneToOneField('protein.ProteinFamily',null=False, on_delete=models.CASCADE, related_name='class_representative_protein_family')
    species = models.ForeignKey('protein.Species', on_delete=models.CASCADE, related_name='class_representative_species')

    objects = models.Manager()  # The default manager.
    custom_objects = CustomClassRepresentativeSpeciesManager()  # The custom manager.

class ClassSimilarityType(IntEnum):
    IDENTITY = 0
    SIMILARITY = 1
    @classmethod
    def choices(cls):
        return [(key.value, key.name) for key in cls]

class ClassSimilarityTie(models.Model):
    class_similarity = models.ForeignKey('alignment.ClassSimilarity',null=False, on_delete=models.CASCADE, related_name='class_similarity_tie_class_similarity')
    protein1 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_tie_protein1')
    protein2 = models.ForeignKey('protein.Protein',null=True, on_delete=models.CASCADE, related_name='class_similarity_tie_protein2')
    type = models.IntegerField(choices=ClassSimilarityType.choices(), default=int(ClassSimilarityType.SIMILARITY))
    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['protein1', 'protein2','type'], name='unique_class_similarity_tie_proteins_type'),   
        ]
    def __str__(self):
        key2value = {}
        for key,value in ClassSimilarityType.choices():
            key2value[key] = value
        return str(self.protein1)+" vs "+str(self.protein2)+": "+key2value[self.type]

