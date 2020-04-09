XCFLAGS=${CFLAGS} \
	-O0 -g -std=c++14 -pedantic \
	-Wall -Wextra -Wfatal-errors -Wno-sequence-point \
	-pipe -fno-omit-frame-pointer -fpermissive \
	-fdata-sections -ffunction-sections -Wl,--gc-sections

XLDFLAGS=${LDFLAGS} \
	-lginac -lcln -ldl -lgmp -lnauty

XCFLAGS_STATIC=${XCFLAGS} -Os -s -static

XLDFLAGS_STATIC=${XLDFLAGS}

all: feynson

feynson: feynson.cpp ginacutils.cpp blake2b.c
	@date "+static const char VERSION[] = \"Feynson $$(hg id -i), built on %Y-%m-%d\n\";" >version.h
	${CXX} ${XCFLAGS} -include version.h -o $@ $< ${XLDFLAGS}

feynson.static: feynson.cpp ginacutils.cpp blake2b.c
	@date "+static const char VERSION[] = \"Feynson $$(hg id -i), built on %Y-%m-%d\n\";" >version.h
	${CXX} ${XCFLAGS_STATIC} -include version.h -o $@ $< ${XLDFLAGS_STATIC}
	@upx --best "$@"

README.md: feynson.cpp mkmanual.sh
	sed '/MANUAL/{n;q}' $@ >$@.tmp
	./mkmanual.sh >>$@.tmp <$<
	mv $@.tmp $@

clean:
	rm -f feynson feynson.static version.h

phony:;
