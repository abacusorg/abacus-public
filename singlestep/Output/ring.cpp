class ringbuffer {
public:
    ringbuffer(int _n) { n = _n; h=0;t=0; s=0; ring = new iorequest[n]; chk();  }

    ~ringbuffer(void) { chk(); delete[] ring; }

    int ringlength(void)  { return s; }
    int isfull(void)  { return (s==n)?1:0; }
    int isempty(void) { return (s==0)?1:0; }
    int isnotempty(void) { return (s==0)?0:1; }

    void pushhead(iorequest r) { 
        assert(!isfull()); h--; if(h<0) h+=n; ring[h] = r; s++; chk();  
    }

    void push(iorequest r) { 
        assert(!isfull());  ring[t] = r; t=(t+1)%n; s++;  chk(); 
    }

    iorequest pop(void) { 
        assert(!isempty()); iorequest r = ring[h]; h=(h+1)%n; s--; chk(); return r; 
    }

    void pr(void) { 
        chk();
        fmt::print("n={:3d}::s={:3d}::h={:3d}::t={:3d}::",n,s,h,t);
        if(s==0) fmt::print("<<empty>>");
        else  {
            for(int i=0;i<s;i++) fmt::print("<{:d}> ", ring[ (h+i)%n ].arenatype);
        }
        fmt::print("\n");
        if(s!=0) {
            for(int i=0;i<s;i++)  ring[ (h+i)%n ].dumpior();
        }
    }

    void chk(void) {
        assert( (h>=0) && (h<n) );
        assert( (t>=0) && (t<n) );
        assert( (s>=0) && (s<=n) );
        if(s==0) assert(h==t);
        assert( (h+s)%n == t ); 
    }
private:
    iorequest *ring;
    int h,t,s,n;
};
