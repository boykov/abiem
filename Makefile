all:
	make --silent -C src testcase tn=.testBIEtest_sigm
	make --silent -C src testcase tn=.testBIEsmall
	make --silent -C src testcase tn=.testBIEsmall3
	make --silent -C src testcase tn=.testBIEsmallNG
	make --silent -C src testcase tn=.testBIEmicro6

test:
	make -C src test

testcase:
	make -C src testcase

clear:
	git clean -xf
