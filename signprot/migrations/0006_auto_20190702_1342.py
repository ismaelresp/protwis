# Generated by Django 2.0.4 on 2019-07-02 11:42

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('signprot', '0005_auto_20180809_1647'),
    ]

    operations = [
        migrations.AlterField(
            model_name='signprotcomplex',
            name='beta_chain',
            field=models.CharField(max_length=1, null=True),
        ),
        migrations.AlterField(
            model_name='signprotcomplex',
            name='beta_protein',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='beta_protein', to='protein.Protein'),
        ),
        migrations.AlterField(
            model_name='signprotcomplex',
            name='gamma_chain',
            field=models.CharField(max_length=1, null=True),
        ),
        migrations.AlterField(
            model_name='signprotcomplex',
            name='gamma_protein',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='gamma_protein', to='protein.Protein'),
        ),
    ]