#!/bin/bash
#
#$ -S /bin/bash
#$ -N trin_eph_indoor
#$ -cwd
#$ -j n
#$ -m e
#$ -M m.brasseur@leibniz-zfmk.de

module load trinity/2.13.2

Trinity --full_cleanup --seqType fq --max_memory 700G --CPU $NSLOTS --output ../stranded_trinity_ephemera_indoor --no_version_check --bflyCalculateCPU --SS_lib_type RF \
--left DE03NGSUKBR124958_1_polyx_trimmed_val_1.fq.gz,DE07NGSUKBR124983_1_polyx_trimmed_val_1.fq.gz,DE08NGSUKBR124965_1_polyx_trimmed_val_1.fq.gz,DE12NGSUKBR124990_1_polyx_trimmed_val_1.fq.gz,DE16NGSUKBR125015_1_polyx_trimmed_val_1.fq.gz,DE17NGSUKBR124997_1_polyx_trimmed_val_1.fq.gz,DE21NGSUKBR125022_1_polyx_trimmed_val_1.fq.gz,DE24NGSUKBR124968_1_polyx_trimmed_val_1.fq.gz,DE25NGSUKBR124950_1_polyx_trimmed_val_1.fq.gz,DE26NGSUKBR125029_1_polyx_trimmed_val_1.fq.gz,DE29NGSUKBR124975_1_polyx_trimmed_val_1.fq.gz,DE30NGSUKBR124957_1_polyx_trimmed_val_1.fq.gz,DE34NGSUKBR124982_1_polyx_trimmed_val_1.fq.gz,DE38NGSUKBR125007_1_polyx_trimmed_val_1.fq.gz,DE39NGSUKBR124989_1_polyx_trimmed_val_1.fq.gz,DE43NGSUKBR125014_1_polyx_trimmed_val_1.fq.gz,DE46NGSUKBR124960_1_polyx_trimmed_val_1.fq.gz,DE48NGSUKBR125021_1_polyx_trimmed_val_1.fq.gz,DE51NGSUKBR124967_1_polyx_trimmed_val_1.fq.gz,DE52NGSUKBR124949_1_polyx_trimmed_val_1.fq.gz,DE56NGSUKBR124974_1_polyx_trimmed_val_1.fq.gz,DE60NGSUKBR124999_1_polyx_trimmed_val_1.fq.gz,DE61NGSUKBR124981_1_polyx_trimmed_val_1.fq.gz,DE65NGSUKBR125006_1_polyx_trimmed_val_1.fq.gz,DE68NGSUKBR124952_1_polyx_trimmed_val_1.fq.gz,DE69NGSUKBR125031_1_polyx_trimmed_val_1.fq.gz,DE70NGSUKBR125013_1_polyx_trimmed_val_1.fq.gz,DE73NGSUKBR124959_1_polyx_trimmed_val_1.fq.gz,DE78NGSUKBR124966_1_polyx_trimmed_val_1.fq.gz,DE82NGSUKBR124991_1_polyx_trimmed_val_1.fq.gz,DE83NGSUKBR124973_1_polyx_trimmed_val_1.fq.gz,DE87NGSUKBR124998_1_polyx_trimmed_val_1.fq.gz,DE91NGSUKBR125023_1_polyx_trimmed_val_1.fq.gz,DE92NGSUKBR125005_1_polyx_trimmed_val_1.fq.gz,DE95NGSUKBR124951_1_polyx_trimmed_val_1.fq.gz,DE96NGSUKBR125030_1_polyx_trimmed_val_1.fq.gz \
--right DE03NGSUKBR124958_2_polyx_trimmed_val_2.fq.gz,DE07NGSUKBR124983_2_polyx_trimmed_val_2.fq.gz,DE08NGSUKBR124965_2_polyx_trimmed_val_2.fq.gz,DE12NGSUKBR124990_2_polyx_trimmed_val_2.fq.gz,DE16NGSUKBR125015_2_polyx_trimmed_val_2.fq.gz,DE17NGSUKBR124997_2_polyx_trimmed_val_2.fq.gz,DE21NGSUKBR125022_2_polyx_trimmed_val_2.fq.gz,DE24NGSUKBR124968_2_polyx_trimmed_val_2.fq.gz,DE25NGSUKBR124950_2_polyx_trimmed_val_2.fq.gz,DE26NGSUKBR125029_2_polyx_trimmed_val_2.fq.gz,DE29NGSUKBR124975_2_polyx_trimmed_val_2.fq.gz,DE30NGSUKBR124957_2_polyx_trimmed_val_2.fq.gz,DE34NGSUKBR124982_2_polyx_trimmed_val_2.fq.gz,DE38NGSUKBR125007_2_polyx_trimmed_val_2.fq.gz,DE39NGSUKBR124989_2_polyx_trimmed_val_2.fq.gz,DE43NGSUKBR125014_2_polyx_trimmed_val_2.fq.gz,DE46NGSUKBR124960_2_polyx_trimmed_val_2.fq.gz,DE48NGSUKBR125021_2_polyx_trimmed_val_2.fq.gz,DE51NGSUKBR124967_2_polyx_trimmed_val_2.fq.gz,DE52NGSUKBR124949_2_polyx_trimmed_val_2.fq.gz,DE56NGSUKBR124974_2_polyx_trimmed_val_2.fq.gz,DE60NGSUKBR124999_2_polyx_trimmed_val_2.fq.gz,DE61NGSUKBR124981_2_polyx_trimmed_val_2.fq.gz,DE65NGSUKBR125006_2_polyx_trimmed_val_2.fq.gz,DE68NGSUKBR124952_2_polyx_trimmed_val_2.fq.gz,DE69NGSUKBR125031_2_polyx_trimmed_val_2.fq.gz,DE70NGSUKBR125013_2_polyx_trimmed_val_2.fq.gz,DE73NGSUKBR124959_2_polyx_trimmed_val_2.fq.gz,DE78NGSUKBR124966_2_polyx_trimmed_val_2.fq.gz,DE82NGSUKBR124991_2_polyx_trimmed_val_2.fq.gz,DE83NGSUKBR124973_2_polyx_trimmed_val_2.fq.gz,DE87NGSUKBR124998_2_polyx_trimmed_val_2.fq.gz,DE91NGSUKBR125023_2_polyx_trimmed_val_2.fq.gz,DE92NGSUKBR125005_2_polyx_trimmed_val_2.fq.gz,DE95NGSUKBR124951_2_polyx_trimmed_val_2.fq.gz,DE96NGSUKBR125030_2_polyx_trimmed_val_2.fq.gz
