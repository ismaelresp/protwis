# Generated by Django 2.0.4 on 2018-05-02 13:50

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
        ('protein', '0002_auto_20180117_1457'),
        ('structure', '0002_structure_sodium'),
    ]

    operations = [
        migrations.CreateModel(
            name='StructureComplexModel',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('pdb', models.TextField()),
                ('version', models.DateField()),
                ('main_template', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='structure.Structure')),
                ('receptor_protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='+', to='protein.Protein')),
                ('sign_protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='+', to='protein.Protein')),
            ],
            options={
                'db_table': 'structure_complex_model',
            },
        ),
        migrations.CreateModel(
            name='StructureComplexModelSeqSim',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('similarity', models.IntegerField()),
                ('homology_model', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='structure.StructureComplexModel')),
                ('template', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='structure.Structure')),
            ],
            options={
                'db_table': 'structure_complex_model_seqsim',
            },
        ),
        migrations.CreateModel(
            name='StructureComplexModelStatsRotamer',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('backbone_template', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='+', to='structure.Structure')),
                ('homology_model', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='structure.StructureComplexModel')),
                ('protein', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='protein.Protein')),
                ('residue', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='residue.Residue')),
                ('rotamer_template', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='+', to='structure.Structure')),
            ],
            options={
                'db_table': 'structure_complex_model_stats_rotamer',
            },
        ),
    ]
