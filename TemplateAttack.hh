#ifndef _TEMPLATE_ATTACK_HEADER
#define _TEMPLATE_ATTACK_HEADER

#define set1 "normal"
#define set2 "default_cache"
#define set3 "cache_associate8"
#define set4 "64kcache"
#define set5 "64k_associate8"
#define set6 "l2default"
#define set7 "l2_l1associate8"
#define set8 "l24m_default"
#define set9 "l24m_l1associate8"
#define set10 "dram_l2default"
#define set11 "dram_l24m_default"
#define set12 "dram_l24m_l1associate8"
#define set13 "dram_l2_l1associate8"


#define _TemplateRange 256
#define _TrainRange 14
#define TemplateRange(i) for(int i = 0; i <_TemplateRange; i++)
#define TrainRange(i) for(int i = 0; i < _TrainRange; i++)
#define TrainSize 28000
#define TestSize 2000
#define Amplify 1000
#define VarianceAmplify 0.5
#define OneRun 3
#define Rounds (_TrainRange + (OneRun - _TrainRange%OneRun)%OneRun)/OneRun // hacky way of ceiling function
#define ColSize 1000
#define elastic_method "average"




#endif
