from django.db import models
from django import forms


# class UploadFileForm(forms.Form):
#     title = forms.CharField(max_length=50)
#     file = forms.FileField()


class unique_id(models.Model):
    UID = models.CharField(primary_key=True, max_length=17)

    def __str__(self):
        return self.UID


class cookie_sequence(models.Model):
    Strand = models.CharField(null=True, max_length=2)
    Spacer = models.CharField(null=True, max_length=50)
    PAM = models.CharField(null=True, max_length=5)
    PBS = models.CharField(null=True, max_length=50)
    RT = models.CharField(null=True, max_length=50)
    EditToNickDistance = models.IntegerField(null=True, blank=False)
    sgRNASpacer = models.CharField(null=True, max_length=50)
    sgRNAPAM = models.CharField(null=True, max_length=50)
    NickToNickDistance = models.IntegerField(null=True, blank=False)
    EditingScore = models.FloatField(null=True, blank=False)

    FID = models.ForeignKey('unique_id', on_delete=models.PROTECT)


class cookie_sequence_index(models.Model):
    UID = models.OneToOneField('unique_id', on_delete=models.PROTECT)
    PAM = models.CharField(max_length=5)
    CUT_SIZE = models.IntegerField(null=True, blank=False)
    MAX_PBS = models.IntegerField(null=True, blank=False)
    MAX_RT = models.IntegerField(null=True, blank=False)
    MIN_PBS = models.IntegerField(null=True, blank=False)
    MIN_RT = models.IntegerField(null=True, blank=False)
    TOP_N = models.IntegerField(null=True, blank=False)
    GENOME = models.CharField(max_length=20)

    MIN_DISGRNA = models.IntegerField(null=True, blank=False)
    MAX_DISGRNA = models.IntegerField(null=True, blank=False)
    HOMOLOGY = models.IntegerField(null=True, blank=False)
    SEQUENCE = models.CharField(max_length=100)


class cookie_position(models.Model):
    Strand = models.CharField(max_length=2)
    Spacer = models.CharField(null=True, max_length=50)
    PAM = models.CharField(null=True, max_length=5)
    PBS = models.CharField(null=True, max_length=50)
    RT = models.CharField(null=True, max_length=50)
    EditToNickDistance = models.IntegerField(null=True, blank=False)
    sgRNASpacer = models.CharField(null=True, max_length=50)
    sgRNAPAM = models.CharField(null=True, max_length=50)
    NickToNickDistance = models.IntegerField(null=True, blank=False)
    EditingScore = models.FloatField(null=True, blank=False)

    FID = models.ForeignKey('unique_id', on_delete=models.PROTECT)


class cookie_position_index(models.Model):
    UID = models.OneToOneField('unique_id', on_delete=models.PROTECT)
    PAM = models.CharField(max_length=5)
    CUT_SIZE = models.IntegerField(null=True, blank=False)
    MAX_PBS = models.IntegerField(null=True, blank=False)
    MAX_RT = models.IntegerField(null=True, blank=False)
    MIN_PBS = models.IntegerField(null=True, blank=False)
    MIN_RT = models.IntegerField(null=True, blank=False)
    TOP_N = models.IntegerField(null=True, blank=False)
    GENOME = models.CharField(max_length=20)

    MIN_DISGRNA = models.IntegerField(null=True, blank=False)
    MAX_DISGRNA = models.IntegerField(null=True, blank=False)
    HOMOLOGY = models.IntegerField(null=True, blank=False)

    Chromosome = models.CharField(null=True, max_length=100)
    Position = models.IntegerField(null=True, blank=False)
    Pattern = models.CharField(null=True, max_length=50)


class cookie_opedvar(models.Model):
    AlleleID = models.IntegerField(null=True, blank=False)
    Type = models.CharField(null=True, max_length=50)
    GeneID = models.CharField(null=True, max_length=50)
    GeneSymbol = models.CharField(null=True, max_length=50)
    HGNC_ID = models.IntegerField(null=True, blank=False)
    Chromosome = models.IntegerField(null=True, blank=False)
    Start = models.IntegerField(null=True, blank=False)
    Stop = models.IntegerField(null=True, blank=False)
    ReferenceAllele = models.CharField(null=True, max_length=5)
    AlternateAllele = models.CharField(null=True, max_length=5)
    
    Strand = models.CharField(null=True, max_length=2)

    Spacer = models.CharField(max_length=50)

    PBS = models.CharField(null=True, max_length=10)
    RT = models.CharField(null=True, max_length=20)
    EditToNickDistance = models.IntegerField(null=True, blank=False)
    sgRNASpacer = models.CharField(null=True, max_length=50)
    NickToNickDistance = models.IntegerField(null=True, blank=False)

    FID = models.ForeignKey('unique_id', on_delete=models.PROTECT)


class cookie_opedvar_index(models.Model):
    UID = models.OneToOneField('unique_id', on_delete=models.PROTECT)
    queryType = models.CharField(max_length=50)
    queryItem = models.CharField(max_length=50)
    PAM = models.CharField(max_length=15)
    GENOME = models.CharField(max_length=20)
    DIRECTION = models.CharField(max_length=15)






