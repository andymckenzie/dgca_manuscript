#need to be in the directory of sim_roc_classes_classes.R for this script to work
#run it via make -f sim_roc_classes_classes.mk
MAKEFLAGS += -j 6  # default number of parallel jobs

.PHONY: all

all: n100f050_1 n100f050_2 n100f050_3 n100f050_4 n100f050_5 \
	n200f050_1 n200f050_2 n200f050_3 n200f050_4 n200f050_5 \
	n300f050_1 n300f050_2 n300f050_3 n300f050_4 n300f050_5 \
	n400f050_1 n400f050_2 n400f050_3 n400f050_4 n400f050_5 \
	n100f060_1 n100f060_2 n100f060_3 n100f060_4 n100f060_5 \
	n200f060_1 n200f060_2 n200f060_3 n200f060_4 n200f060_5 \
	n300f060_1 n300f060_2 n300f060_3 n300f060_4 n300f060_5 \
	n400f060_1 n400f060_2 n400f060_3 n400f060_4 n400f060_5 \
	n100f070_1 n100f070_2 n100f070_3 n100f070_4 n100f070_5 \
	n200f070_1 n200f070_2 n200f070_3 n200f070_4 n200f070_5 \
	n300f070_1 n300f070_2 n300f070_3 n300f070_4 n300f070_5 \
	n400f070_1 n400f070_2 n400f070_3 n400f070_4 n400f070_5 \
	n100f080_1 n100f080_2 n100f080_3 n100f080_4 n100f080_5 \
	n200f080_1 n200f080_2 n200f080_3 n200f080_4 n200f080_5 \
	n300f080_1 n300f080_2 n300f080_3 n300f080_4 n300f080_5 \
	n400f080_1 n400f080_2 n400f080_3 n400f080_4 n400f080_5 \
	n100f090_1 n100f090_2 n100f090_3 n100f090_4 n100f090_5 \
	n200f090_1 n200f090_2 n200f090_3 n200f090_4 n200f090_5 \
	n300f090_1 n300f090_2 n300f090_3 n300f090_4 n300f090_5 \
	n400f090_1 n400f090_2 n400f090_3 n400f090_4 n400f090_5 \
	n100f100_1 n100f100_2 n100f100_3 n100f100_4 n100f100_5 \
	n200f100_1 n200f100_2 n200f100_3 n200f100_4 n200f100_5 \
	n300f100_1 n300f100_2 n300f100_3 n300f100_4 n300f100_5 \
	n400f100_1 n400f100_2 n400f100_3 n400f100_4 n400f100_5

######

n100f050_1: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.50" "TRUE" "TRUE" "TRUE"
n100f050_2: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.50" "TRUE" "TRUE" "TRUE"
n100f050_3: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.50" "TRUE" "TRUE" "TRUE"
n100f050_4: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.50" "TRUE" "TRUE" "TRUE"
n100f050_5: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.50" "TRUE" "TRUE" "TRUE"

n200f050_1: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.50" "TRUE" "TRUE" "TRUE"
n200f050_2: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.50" "TRUE" "TRUE" "TRUE"
n200f050_3: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.50" "TRUE" "TRUE" "TRUE"
n200f050_4: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.50" "TRUE" "TRUE" "TRUE"
n200f050_5: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.50" "TRUE" "TRUE" "TRUE"

n300f050_1: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.50" "TRUE" "TRUE" "TRUE"
n300f050_2: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.50" "TRUE" "TRUE" "TRUE"
n300f050_3: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.50" "TRUE" "TRUE" "TRUE"
n300f050_4: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.50" "TRUE" "TRUE" "TRUE"
n300f050_5: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.50" "TRUE" "TRUE" "TRUE"

n400f050_1: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.50" "TRUE" "TRUE" "TRUE"
n400f050_2: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.50" "TRUE" "TRUE" "TRUE"
n400f050_3: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.50" "TRUE" "TRUE" "TRUE"
n400f050_4: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.50" "TRUE" "TRUE" "TRUE"
n400f050_5: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.50" "TRUE" "TRUE" "TRUE"

#####

n100f060_1: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.60" "TRUE" "TRUE" "TRUE"
n100f060_2: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.60" "TRUE" "TRUE" "TRUE"
n100f060_3: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.60" "TRUE" "TRUE" "TRUE"
n100f060_4: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.60" "TRUE" "TRUE" "TRUE"
n100f060_5: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.60" "TRUE" "TRUE" "TRUE"

n200f060_1: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.60" "TRUE" "TRUE" "TRUE"
n200f060_2: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.60" "TRUE" "TRUE" "TRUE"
n200f060_3: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.60" "TRUE" "TRUE" "TRUE"
n200f060_4: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.60" "TRUE" "TRUE" "TRUE"
n200f060_5: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.60" "TRUE" "TRUE" "TRUE"

n300f060_1: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.60" "TRUE" "TRUE" "TRUE"
n300f060_2: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.60" "TRUE" "TRUE" "TRUE"
n300f060_3: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.60" "TRUE" "TRUE" "TRUE"
n300f060_4: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.60" "TRUE" "TRUE" "TRUE"
n300f060_5: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.60" "TRUE" "TRUE" "TRUE"

n400f060_1: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.60" "TRUE" "TRUE" "TRUE"
n400f060_2: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.60" "TRUE" "TRUE" "TRUE"
n400f060_3: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.60" "TRUE" "TRUE" "TRUE"
n400f060_4: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.60" "TRUE" "TRUE" "TRUE"
n400f060_5: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.60" "TRUE" "TRUE" "TRUE"

#####

n100f070_1: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.70" "TRUE" "TRUE" "TRUE"
n100f070_2: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.70" "TRUE" "TRUE" "TRUE"
n100f070_3: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.70" "TRUE" "TRUE" "TRUE"
n100f070_4: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.70" "TRUE" "TRUE" "TRUE"
n100f070_5: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.70" "TRUE" "TRUE" "TRUE"

n200f070_1: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.70" "TRUE" "TRUE" "TRUE"
n200f070_2: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.70" "TRUE" "TRUE" "TRUE"
n200f070_3: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.70" "TRUE" "TRUE" "TRUE"
n200f070_4: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.70" "TRUE" "TRUE" "TRUE"
n200f070_5: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.70" "TRUE" "TRUE" "TRUE"

n300f070_1: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.70" "TRUE" "TRUE" "TRUE"
n300f070_2: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.70" "TRUE" "TRUE" "TRUE"
n300f070_3: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.70" "TRUE" "TRUE" "TRUE"
n300f070_4: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.70" "TRUE" "TRUE" "TRUE"
n300f070_5: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.70" "TRUE" "TRUE" "TRUE"

n400f070_1: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.70" "TRUE" "TRUE" "TRUE"
n400f070_2: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.70" "TRUE" "TRUE" "TRUE"
n400f070_3: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.70" "TRUE" "TRUE" "TRUE"
n400f070_4: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.70" "TRUE" "TRUE" "TRUE"
n400f070_5: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.70" "TRUE" "TRUE" "TRUE"

#####

n100f080_1: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.80" "TRUE" "TRUE" "TRUE"
n100f080_2: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.80" "TRUE" "TRUE" "TRUE"
n100f080_3: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.80" "TRUE" "TRUE" "TRUE"
n100f080_4: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.80" "TRUE" "TRUE" "TRUE"
n100f080_5: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.80" "TRUE" "TRUE" "TRUE"

n200f080_1: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.80" "TRUE" "TRUE" "TRUE"
n200f080_2: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.80" "TRUE" "TRUE" "TRUE"
n200f080_3: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.80" "TRUE" "TRUE" "TRUE"
n200f080_4: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.80" "TRUE" "TRUE" "TRUE"
n200f080_5: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.80" "TRUE" "TRUE" "TRUE"

n300f080_1: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.80" "TRUE" "TRUE" "TRUE"
n300f080_2: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.80" "TRUE" "TRUE" "TRUE"
n300f080_3: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.80" "TRUE" "TRUE" "TRUE"
n300f080_4: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.80" "TRUE" "TRUE" "TRUE"
n300f080_5: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.80" "TRUE" "TRUE" "TRUE"

n400f080_1: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.80" "TRUE" "TRUE" "TRUE"
n400f080_2: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.80" "TRUE" "TRUE" "TRUE"
n400f080_3: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.80" "TRUE" "TRUE" "TRUE"
n400f080_4: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.80" "TRUE" "TRUE" "TRUE"
n400f080_5: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.80" "TRUE" "TRUE" "TRUE"

#####

n100f090_1: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.90" "TRUE" "TRUE" "TRUE"
n100f090_2: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.90" "TRUE" "TRUE" "TRUE"
n100f090_3: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.90" "TRUE" "TRUE" "TRUE"
n100f090_4: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.90" "TRUE" "TRUE" "TRUE"
n100f090_5: sim_modules.R
	Rscript sim_modules.R "100" "random" "0.90" "TRUE" "TRUE" "TRUE"

n200f090_1: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.90" "TRUE" "TRUE" "TRUE"
n200f090_2: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.90" "TRUE" "TRUE" "TRUE"
n200f090_3: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.90" "TRUE" "TRUE" "TRUE"
n200f090_4: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.90" "TRUE" "TRUE" "TRUE"
n200f090_5: sim_modules.R
	Rscript sim_modules.R "200" "random" "0.90" "TRUE" "TRUE" "TRUE"

n300f090_1: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.90" "TRUE" "TRUE" "TRUE"
n300f090_2: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.90" "TRUE" "TRUE" "TRUE"
n300f090_3: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.90" "TRUE" "TRUE" "TRUE"
n300f090_4: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.90" "TRUE" "TRUE" "TRUE"
n300f090_5: sim_modules.R
	Rscript sim_modules.R "300" "random" "0.90" "TRUE" "TRUE" "TRUE"

n400f090_1: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.90" "TRUE" "TRUE" "TRUE"
n400f090_2: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.90" "TRUE" "TRUE" "TRUE"
n400f090_3: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.90" "TRUE" "TRUE" "TRUE"
n400f090_4: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.90" "TRUE" "TRUE" "TRUE"
n400f090_5: sim_modules.R
	Rscript sim_modules.R "400" "random" "0.90" "TRUE" "TRUE" "TRUE"

#####

n100f100_1: sim_modules.R
	Rscript sim_modules.R "100" "random" "1.00" "TRUE" "TRUE" "TRUE"
n100f100_2: sim_modules.R
	Rscript sim_modules.R "100" "random" "1.00" "TRUE" "TRUE" "TRUE"
n100f100_3: sim_modules.R
	Rscript sim_modules.R "100" "random" "1.00" "TRUE" "TRUE" "TRUE"
n100f100_4: sim_modules.R
	Rscript sim_modules.R "100" "random" "1.00" "TRUE" "TRUE" "TRUE"
n100f100_5: sim_modules.R
	Rscript sim_modules.R "100" "random" "1.00" "TRUE" "TRUE" "TRUE"

n200f100_1: sim_modules.R
	Rscript sim_modules.R "200" "random" "1.00" "TRUE" "TRUE" "TRUE"
n200f100_2: sim_modules.R
	Rscript sim_modules.R "200" "random" "1.00" "TRUE" "TRUE" "TRUE"
n200f100_3: sim_modules.R
	Rscript sim_modules.R "200" "random" "1.00" "TRUE" "TRUE" "TRUE"
n200f100_4: sim_modules.R
	Rscript sim_modules.R "200" "random" "1.00" "TRUE" "TRUE" "TRUE"
n200f100_5: sim_modules.R
	Rscript sim_modules.R "200" "random" "1.00" "TRUE" "TRUE" "TRUE"

n300f100_1: sim_modules.R
	Rscript sim_modules.R "300" "random" "1.00" "TRUE" "TRUE" "TRUE"
n300f100_2: sim_modules.R
	Rscript sim_modules.R "300" "random" "1.00" "TRUE" "TRUE" "TRUE"
n300f100_3: sim_modules.R
	Rscript sim_modules.R "300" "random" "1.00" "TRUE" "TRUE" "TRUE"
n300f100_4: sim_modules.R
	Rscript sim_modules.R "300" "random" "1.00" "TRUE" "TRUE" "TRUE"
n300f100_5: sim_modules.R
	Rscript sim_modules.R "300" "random" "1.00" "TRUE" "TRUE" "TRUE"

n400f100_1: sim_modules.R
	Rscript sim_modules.R "400" "random" "1.00" "TRUE" "TRUE" "TRUE"
n400f100_2: sim_modules.R
	Rscript sim_modules.R "400" "random" "1.00" "TRUE" "TRUE" "TRUE"
n400f100_3: sim_modules.R
	Rscript sim_modules.R "400" "random" "1.00" "TRUE" "TRUE" "TRUE"
n400f100_4: sim_modules.R
	Rscript sim_modules.R "400" "random" "1.00" "TRUE" "TRUE" "TRUE"
n400f100_5: sim_modules.R
	Rscript sim_modules.R "400" "random" "1.00" "TRUE" "TRUE" "TRUE"
