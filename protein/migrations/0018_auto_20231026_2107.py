# Generated by Django 3.0.3 on 2023-10-26 19:07

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0017_auto_20231025_0044'),
    ]

    operations = [
        migrations.AlterField(
            model_name='proteincouplings',
            name='deltaGDP_conc_family',
            field=models.DecimalField(decimal_places=2, max_digits=4, null=True),
        ),
    ]
