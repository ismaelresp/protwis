from django.db import models
from django.db import connection
# Create your models here.

class AlignmentConsensus(models.Model):
    slug = models.SlugField(max_length=100, unique=True)
    alignment = models.BinaryField()
    gn_consensus = models.BinaryField(blank=True) # Store conservation calculation for each GN


class CustomClassSimilarityManager(models.Manager):
    def truncate_table(self):
        cursor = connection.cursor()
        table_name = self.model._meta.db_table
        sql = 'TRUNCATE TABLE "{0}"'.format(table_name)
        cursor.execute(sql)

class ClassSimilarity(models.Model):
    protein_family1 = models.ForeignKey('protein.ProteinFamily',null=False, on_delete=models.CASCADE, related_name='class_similarity_protein_family1')
    protein_family2 = models.ForeignKey('protein.ProteinFamily',null=False, on_delete=models.CASCADE, related_name='class_similarity_protein_family2')
    protein1 = models.ForeignKey('protein.Protein',null=False, on_delete=models.CASCADE, related_name='class_similarity_protein1')
    protein2 = models.ForeignKey('protein.Protein',null=False, on_delete=models.CASCADE, related_name='class_similarity_protein2')
    similarity = models.IntegerField(null=False)

    objects = models.Manager()  # The default manager.
    custom_objects = CustomClassSimilarityManager()  # The custom manager.
    class Meta:
        constraints = [
            models.UniqueConstraint(fields=['protein_family1', 'protein_family2'], name='unique_class_similarity_protein_family'),
            models.UniqueConstraint(fields=['protein1', 'protein2'], name='unique_class_similarity_protein'),
        ]
    def __str__(self):
        return str(self.protein_family1)+" - "+str(self.protein_family2)+": "+str(self.similarity)
