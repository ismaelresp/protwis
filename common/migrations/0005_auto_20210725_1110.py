# Generated by Django 3.0.3 on 2021-07-25 09:10

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('common', '0004_auto_20210725_1015'),
    ]

    operations = [
        migrations.AlterField(
            model_name='citation',
            name='main',
            field=models.TextField(null=True),
        ),
    ]