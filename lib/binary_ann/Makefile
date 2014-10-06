CXX = /opt/local/bin/g++-mp-4.4 
BIN= "./bin"
CXXFLAGS = -O2 -Wall -fmessage-length=0 -ftree-vectorize -mpopcnt -I./src/
SOURCES=testANN.cpp
TARGET=testANN

build: $(TARGET)

debug: CXX += -DDEBUG -g
debug: $(TARGET)

clean:
	rm -f $(BIN)/$(TARGET)
	
rebuild: clean build
	
$(TARGET):
	$(CXX) $(CXXFLAGS) $(TARGET).cpp -o $(TARGET)
	mkdir -p $(BIN)
	mv $(TARGET) $(BIN)/$(TARGET)
	