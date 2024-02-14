from django.db import models

# Create your models here.

class AlignmentConsensus(models.Model):
    slug = models.SlugField(max_length=100, unique=True)
    alignment = models.BinaryField()
    gn_consensus = models.BinaryField(blank=True) # Store conservation calculation for each GN

class ClassSimilarity(models.Model):
    ProteinFamily1 = models.ForeignKey('protein.ProteinFamily', on_delete=models.CASCADE, related_name='class_similarity_protein_family1')
    ProteinFamily2 = models.ForeignKey('protein.ProteinFamily', on_delete=models.CASCADE, related_name='class_similarity_protein_family2')
    Similarity = models.IntegerField(null=False)
    Protein1 = models.ForeignKey('protein.Protein', on_delete=models.CASCADE, related_name='class_similarity_protein1')
    Protein2 = models.ForeignKey('protein.Protein', on_delete=models.CASCADE, related_name='class_similarity_protein2')
