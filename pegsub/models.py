from django.db import models

# Create your models here.
# 定义预测提交任务的模型类


# 单个序列求pegRNA
# class preInfo(models.Model):
#     Target = models.CharField(max_length=200)
#     PBS = models.CharField(max_length=50)
#     RT = models.CharField(max_length=50)
#     result = models.CharField(max_length=50)
#
#
# # 文件求pegRNA
# class preFileInfo(models.Model):
#     File = models.FileField(null=True, blank=True)
#     result = models.IntegerField(null=True, blank=False)
#
#
# # research 师兄那个求字段的功能
# # PegRNAID	Target(47bp)	Spacer	PBS	RT	NickingCoordinate
# class pegRNAInfo(models.Model):
#     PegRNAID = models.IntegerField(null=True, blank=False)
#     Target = models.CharField(max_length=50)
#     Spacer = models.CharField(max_length=200)
#     PBS = models.CharField(max_length=200)
#     RT = models.CharField(max_length=200)
#     NickingCoordinate = models.CharField(max_length=200)
#
#
# # database  search
#
# class Final_correct_info(models.Model):
#     AlleleID = models.IntegerField(null=True, blank=False)
#     Type = models.CharField(max_length=50)
#     Chromosome = models.CharField(max_length=5)
#     Start = models.IntegerField(null=True, blank=False)
#     Stop = models.IntegerField(null=True, blank=False)
#     ReferenceAllele = models.CharField(max_length=5)
#     AlternateAllele = models.CharField(max_length=5)
#     ReferenceSequence = models.CharField(max_length=200)
#     PegRNAID = models.IntegerField(null=True, blank=False)
#     Target = models.CharField(max_length=50)
#     Spacer = models.CharField(max_length=200)
#     PBS = models.CharField(max_length=50)
#     RT = models.CharField(max_length=50)
#     NickingCoordinate = models.IntegerField(null=True, blank=False)
#     Efficiency = models.FloatField(null=True, blank=False)
#     PAM = models.CharField(max_length=5)
#
#
# class Final_install_info(models.Model):
#     AlleleID = models.IntegerField(null=True, blank=False)
#     Type = models.CharField(max_length=50)
#     Chromosome = models.CharField(max_length=5)
#     Start = models.IntegerField(null=True, blank=False)
#     Stop = models.IntegerField(null=True, blank=False)
#     ReferenceAllele = models.CharField(max_length=5)
#     AlternateAllele = models.CharField(max_length=5)
#     ReferenceSequence = models.CharField(max_length=200)
#     PegRNAID = models.IntegerField(null=True, blank=False)
#     Target = models.CharField(max_length=50)
#     Spacer = models.CharField(max_length=50)
#     PBS = models.CharField(max_length=50)
#     RT = models.CharField(max_length=50)
#     NickingCoordinate = models.IntegerField(null=True, blank=False)
#     Efficiency = models.FloatField(null=True, blank=False)
#     PAM = models.CharField(max_length=5)


"""
# 存储sequence的历史模块 两个模型一个索引模型 一个数据表模型
#展示的结果字段
"""


class cookie_sequence(models.Model):
    # PegRNAID = models.IntegerField(null=True, blank=False)
    Target = models.CharField(max_length=50)
    # Spacer = models.CharField(max_length=50)
    PBS = models.CharField(max_length=50)
    RT = models.CharField(max_length=50)
    # NickingCoordinate = models.IntegerField(null=True, blank=False)
    Efficiency = models.FloatField(null=True, blank=False)
    # Spacer = models.CharField(max_length=50)
    time = models.CharField(max_length=20)



# 索引模型包含sequence模块的参数字段
class cookie_sequence_index(models.Model):
    PAM = models.CharField(max_length=5)
    CUT_SIZE = models.IntegerField(null=True, blank=False)
    MAX_PBS = models.IntegerField(null=True, blank=False)
    MAX_RT = models.IntegerField(null=True, blank=False)
    MIN_PBS = models.IntegerField(null=True, blank=False)
    MIN_RT = models.IntegerField(null=True, blank=False)
    number_show = models.IntegerField(null=True, blank=False)
    input_sequence = models.CharField(max_length=100)
    # index_time = models.ForeignKey(cookie_sequence, models.CASCADE)
    time = models.CharField(max_length=20)


'''
sequence加一个四月的输出内容模型和参数模型 后缀为April 重新注释模型类
'''
#结果模型
class cookie_sequence_April(models.Model):
    Strand=models.CharField(max_length=2)
    Spacer=models.CharField(max_length=50)
    PAM = models.CharField(max_length=5)
    PBS=models.CharField(max_length=50)
    RT=models.CharField(max_length=50)
    EditToNickDistance=models.IntegerField(null=True, blank=False)
    sgRNASpacer=models.CharField(max_length=50)
    sgRNAPAM=models.CharField(max_length=50)
    NickToNickDistance = models.IntegerField(null=True, blank=False)
    Efficiency = models.FloatField(null=True, blank=False)
    time = models.CharField(max_length=20)
###################
#参数模型
class cookie_sequence_April_index(models.Model):
    PAM = models.CharField(max_length=5)
    CUT_SIZE = models.IntegerField(null=True, blank=False)
    MAX_PBS = models.IntegerField(null=True, blank=False)
    MAX_RT = models.IntegerField(null=True, blank=False)
    MIN_PBS = models.IntegerField(null=True, blank=False)
    MIN_RT = models.IntegerField(null=True, blank=False)
    number_show = models.IntegerField(null=True, blank=False)
    ###
    MIN_DISGRNA=models.IntegerField(null=True, blank=False)
    MAX_DISGRNA = models.IntegerField(null=True, blank=False)
    HOMOLOGY= models.IntegerField(null=True, blank=False)
    input_sequence = models.CharField(max_length=100)
    # index_time = models.ForeignKey(cookie_sequence, models.CASCADE)
    time = models.CharField(max_length=20)

"""
position模块的两个模型
cookie_position 数据表模型
cookie_position_index 数据表索引模型
"""


class cookie_position(models.Model):
    # PegRNAID = models.IntegerField(null=True, blank=False)
    Target = models.CharField(max_length=50)
    # Spacer = models.CharField(max_length=50)
    PBS = models.CharField(max_length=50)
    RT = models.CharField(max_length=50)
    # NickingCoordinate = models.IntegerField(null=True, blank=False)
    Efficiency = models.FloatField(null=True, blank=False)
    # Spacer = models.CharField(max_length=50)
    time = models.CharField(max_length=20)


class cookie_position_index(models.Model):
    PAM = models.CharField(max_length=5)
    CUT_SIZE = models.IntegerField(null=True, blank=False)
    MAX_PBS = models.IntegerField(null=True, blank=False)
    MAX_RT = models.IntegerField(null=True, blank=False)
    MIN_PBS = models.IntegerField(null=True, blank=False)
    MIN_RT = models.IntegerField(null=True, blank=False)
    number_show = models.IntegerField(null=True, blank=False)
    input_chr = models.CharField(max_length=100)
    input_position = models.IntegerField(null=True, blank=False)
    input_pattern = models.CharField(max_length=50)
    time = models.CharField(max_length=20)
    # index_time = models.ForeignKey(cookie_position, models.CASCADE)
'''
Position加一个四月的输出内容模型和参数模型 后缀为April 重新注释模型类
'''
class cookie_position_April(models.Model):
    Strand=models.CharField(max_length=2)
    Spacer=models.CharField(max_length=50)
    PAM = models.CharField(max_length=5)
    PBS=models.CharField(max_length=50)
    RT=models.CharField(max_length=50)
    EditToNickDistance=models.IntegerField(null=True, blank=False)
    sgRNASpacer=models.CharField(max_length=50)
    sgRNAPAM=models.CharField(max_length=50)
    NickToNickDistance = models.IntegerField(null=True, blank=False)
    Efficiency = models.FloatField(null=True, blank=False)
    time = models.CharField(max_length=20)

class cookie_position_April_index(models.Model):
    PAM = models.CharField(max_length=5)
    CUT_SIZE = models.IntegerField(null=True, blank=False)
    MAX_PBS = models.IntegerField(null=True, blank=False)
    MAX_RT = models.IntegerField(null=True, blank=False)
    MIN_PBS = models.IntegerField(null=True, blank=False)
    MIN_RT = models.IntegerField(null=True, blank=False)
    number_show = models.IntegerField(null=True, blank=False)
    ###四月更新字段
    MIN_DISGRNA=models.IntegerField(null=True, blank=False)
    MAX_DISGRNA = models.IntegerField(null=True, blank=False)
    HOMOLOGY= models.IntegerField(null=True, blank=False)
    ###四月更新字段end
    input_chr = models.CharField(max_length=100)
    input_position = models.IntegerField(null=True, blank=False)
    input_pattern = models.CharField(max_length=50)
    time = models.CharField(max_length=20)





"""
data_search 模块的两个模型
cookie_position 数据表模型
cookie_position_index 数据表索引模型
"""


class cookie_databse(models.Model):
    AlleleID = models.IntegerField(null=True, blank=False)
    Type = models.CharField(max_length=50)
    Chromosome = models.IntegerField(null=True, blank=False)
    Start = models.IntegerField(null=True, blank=False)
    Stop = models.IntegerField(null=True, blank=False)
    ReferenceAllele = models.CharField(max_length=5)
    AlternateAllele = models.CharField(max_length=5)
    Spacer = models.CharField(max_length=50)
    PBS = models.CharField(max_length=10)
    RT = models.CharField(max_length=20)
    # PAM_choice = (
    #     (0, 'NG'),
    #     (1, 'NGG'),
    # )
    # PAM = models.IntegerField(choices=PAM_choice, default=0)
    PAM = models.CharField(max_length=5)
    # TIME = models.DateTimeField()
    time = models.CharField(max_length=10)


class cookie_database_index(models.Model):
    ALLELEID = models.IntegerField(null=True, blank=False)
    # PAM_choice = (
    #     (0, 'NG'),
    #     (1, 'NGG'),
    # )
    # PAM = models.IntegerField(choices=PAM_choice, default=0)
    PAM = models.CharField(max_length=5)
    # DIRECTION_choice = (
    #     (0, 'install'),
    #     (1, 'direction')
    # )
    # DIRECTION = models.IntegerField(choices=DIRECTION_choice, default=1)
    DIRECTION = models.CharField(max_length=15)
    # TIME = models.DateTimeField()
    time = models.CharField(max_length=10)
    # index_time = models.ForeignKey(cookie_databse, models.CASCADE)



"""
data_search 模块的两个模型
cookie_position_April 数据表模型
cookie_position_April_index 数据表索引模型
"""
class cookie_database_April(models.Model):
    AlleleID = models.IntegerField(null=True, blank=False)
    Type = models.CharField(max_length=50)
    GeneID=models.CharField(max_length=50)
    GeneSymbol=models.CharField(max_length=50)
    HGNC_ID=models.IntegerField(null=True,blank=False)
    Chromosome = models.IntegerField(null=True, blank=False)
    Start = models.IntegerField(null=True, blank=False)
    Stop = models.IntegerField(null=True, blank=False)
    ReferenceAllele = models.CharField(max_length=5)
    AlternateAllele = models.CharField(max_length=5)
    Spacer = models.CharField(max_length=50)
    PBS = models.CharField(max_length=10)
    RT = models.CharField(max_length=20)
    EditToNickDistance=models.IntegerField(null=True,blank=False)

    sgRNASpacer=models.CharField(max_length=50)

    NickToNickDistance=models.IntegerField(null=True,blank=False)

    # TIME = models.DateTimeField()
    PAM=models.CharField(max_length=15)
    #
    time = models.CharField(max_length=10)

class cookie_database_April_index(models.Model):
    queryType=models.CharField(max_length=50)
    queryItem=models.IntegerField(null=True,blank=False)
    PAM=models.CharField(max_length=15)
    DIRECTION = models.CharField(max_length=15)
    time = models.CharField(max_length=10)