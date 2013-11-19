# ======================
# = top level Makefile =
# ======================
DIRS = ProgressBar TruthNtupleLooper mt2 TruthRecordHelpers BMinusLCutflow \
	   HistogramHandlers

all::build link
# all::build_and_link

build_and_link::
	@ for dir in $(DIRS); \
	do (cd $$dir ; echo "" ; echo "Building $$dir" ; echo "" ; make shlib ; echo "" ; echo "Linking $$dir" ; echo "" ; make executable ); \
	done

build::
	@ for dir in $(DIRS); \
	do (cd $$dir ; echo "" ; echo "Building $$dir" ; echo "" ; make shlib ); \
	done
link::
	@ for dir in $(DIRS); \
	do (cd $$dir ; echo "" ; echo "Linking $$dir" ; echo "" ; make executable ); \
	done

clean::
	@ for dir in $(DIRS); \
	do (cd $$dir ; make clean ); \
	done
