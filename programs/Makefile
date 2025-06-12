TARGET      := go
CPP         := g++
CC          := gcc
CFLAGS      := -w -MMD -std=c++11 
LIBS        := -lm -lpthread -ltcmalloc -lstdc++
LIBPATH     := ../engine
LDLIBS      := #

PROTOC      := protoc
PROTO_FLAGS := #
PROTO_LIBS  := -lprotobuf 

RM          := rm -rf
OBJS        := recordtype.pb.o $(LIBPATH)/util.o kparts.o

OPTS        := #from command line (eg: OPTS+=...)

-include    $(OBJS:%.o=%.d)

all: CFLAGS += -O3
all: $(TARGET) 

debug: CFLAGS += -DDEBUG -g3
debug: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OPTS) $^ $(LDLIBS) $(LIBS) $(PROTO_LIBS) $(INC) -o $@


%.pb.cc: %.proto
	$(PROTOC) -I=. --cpp_out=. $<

%.pb.o: %.pb.cc %.pb.h
	$(CC) $(CFLAGS) $(OPTS) -c $< -o $@

%.o: %.cpp %.h
	$(CC) $(CFLAGS) $(OPTS) -c $< -o $@

%.o: %.cpp
	$(CC) $(CFLAGS) $(OPTS) -c $< -o $@


clean:
	$(RM) $(LIBPATH)/*.o *.o $(LIBPATH)/*.d *.d $(TARGET) *.pb.*
#	$(RM) $(LIBPATH)/*.o *.o $(LIBPATH)/*.d *.d $(TARGET) $(TARGET_ONE) *.pb.*

# Dependencies:
# 1. Protocol buffers v2.5.0
# 2. Google's tcmalloc v2.1 --> lubunwind 0.99-beta

