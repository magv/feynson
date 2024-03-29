XCFLAGS=${CFLAGS} \
	-O3 -g -std=c++14 -pedantic \
	-Wall -Wextra -Wfatal-errors -Wno-sequence-point \
	-pipe -fno-omit-frame-pointer -fpermissive

XLDFLAGS=${LDFLAGS} \
	-lginac -lcln -ldl -lgmp -lnauty

XCFLAGS_STATIC=${XCFLAGS} -Os -s -static \
	-fdata-sections -ffunction-sections -Wl,--gc-sections

XLDFLAGS_STATIC=${XLDFLAGS}

MATH?=math

all: feynson

feynson: feynson.cpp ginacutils.cpp blake2b.c
	@date "+static const char VERSION[] = \"Feynson $$(git --git-dir=.git rev-parse --short=12 HEAD 2>/dev/null || hg -R. id -i), built on %Y-%m-%d\n\";" >version.h
	${CXX} ${XCFLAGS} -include version.h -o $@ $< ${XLDFLAGS}

feynson.static: feynson.cpp ginacutils.cpp blake2b.c
	env CXX="${CXX}" ./mkversion.sh >version.h
	${CXX} ${XCFLAGS_STATIC} -include version.h -o $@ $< ${XLDFLAGS_STATIC}
	@upx --best "$@"

README.md: feynson.cpp mkmanual.sh
	sed '/MANUAL/{n;q}' $@ >$@.tmp
	./mkmanual.sh >>$@.tmp <$<
	mv $@.tmp $@

test: feynson phony
	${MATH} -script test.m

clean: phony
	rm -f feynson feynson.static version.h

phony:;
