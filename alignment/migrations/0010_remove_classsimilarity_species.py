# Generated by Django 2.2.17 on 2024-03-16 01:55

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('alignment', '0009_classsimilarity_species'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='classsimilarity',
            name='species',
        ),
    ]
