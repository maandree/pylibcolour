.POSIX:

all: lib.py

lib.py: pylibcolour.py matrices.py
	sed '/^### MATRICES %%%$$/q' < pylibcolour.py | sed '$$d' > lib.py
	./matrices.py >> lib.py
	sed '1,/^### MATRICES %%%$$/d' < pylibcolour.py >> lib.py

check: lib.py
	./test.py

clean:
	-rm lib.py

.PHONY: all check clean
