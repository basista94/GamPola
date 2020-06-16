#ifndef GSIMPLERANDOMENGINE_HH
#define GSIMPLERANDOMENGINE_HH

class GSimpleRandomEngine {

public:

    GSimpleRandomEngine(){
	_next=1;
    }

    void reset() {
	_next=1;
    }

    virtual double random();

private:

  unsigned long int _next;

};

#endif


