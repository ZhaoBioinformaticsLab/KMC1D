SOURCE = KMC1D.cu
OBJECT = KMC1D.o
TARGET = KMC1D
CC = nvcc#g++

CMP_OPT  = -O2 
LINK_OPT = 

$(TARGET):$(OBJECT)
	$(CC) $(LINK_OPT) -o $(TARGET) $(OBJECT)

$(OBJECT):$(SOURCE)
	$(CC) $(CMP_OPT) -c $(SOURCE)

clean:
	rm -f $(TARGET) $(OBJECT) 
