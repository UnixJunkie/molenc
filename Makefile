.PHONY: build install uninstall reinstall test

build:
	dune build @install -j `getconf _NPROCESSORS_ONLN`

clean:
	rm -rf _build

edit:
	emacs src/*.ml TODO commands.sh &

install: build
	dune install

uninstall:
	dune uninstall

reinstall: uninstall install

test:
	rm -f _build/default/src/fp_test.exe
	dune build src/fp_test.exe
	_build/default/src/fp_test.exe
