SHELL=/bin/bash
.PHONY: help clone build shell
.DEFAULT_GOAL: help

help:
	@echo "Commands: make link"

build:
	docker build --rm -t dbgwas .

link:
	docker run -v /tmp/data:/Output -it dbgwas "bash"

copy:
	docker-machine scp -r $${machineName}:/tmp/data/ Output/

amikacin:
	# docker run --detach -v /tmp/data:/Output dbgwas bash -c 'bash /Analysis/Amikacin/CALDERA.sh > Output/log 2>&1'
	# docker run --detach -v /tmp/data:/Output dbgwas bash -c 'bash /Analysis/Amikacin/Comp_COIN.sh >> Output/log 2>&1'
	docker run --detach -v /tmp/data:/Output dbgwas bash -c 'bash /Analysis/Amikacin/Comp_AllUnitigs.sh >> Output/log 2>&1'

amuc:
	docker run --detach -v /tmp/data:/Output dbgwas bash -c 'bash /Analysis/Amuc/Caldera.sh > Output/log 2>&1'

ath:
	docker run --detach -v /tmp/data:/Output dbgwas bash -c 'bash /Analysis/Ath/ath.sh > Output/log 2>&1'

explo:
	docker run --detach -v /tmp/data:/Output dbgwas bash -c 'bash /Analysis/Explo/comparisons.sh > Output/log 2>&1'

simu:
	docker run --detach -v /tmp/data:/Output dbgwas bash -c 'bash /Analysis/Simulations/Run_speed.sh >> Output/log 2>&1'
	# docker run --detach -v /tmp/data:/Output dbgwas bash -c 'bash /Analysis/Simulations/Run_imbalance.sh >> Output/log 2>&1'

speed:
	docker run --detach -v /tmp/data:/Output dbgwas bash -c 'bash /Analysis/Simulations/Run_speed_2.sh > Output/log_speed 2>&1'


example:
	'./Analysis/Ath/example.sh > Output/log 2>&1'
