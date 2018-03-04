
#define _CRT_SECURE_NO_DEPRECATE // to keep the VS compiler happy with TBB

// let's use a large type to store fib numbers
typedef unsigned long long fib_type;

#include "tree.h"

// the actual step code computing the fib numbers goes here
int tree_step::execute( const int & tag, tree_context & ctxt ) const
{
    switch( tag ) {
        case 0 : ctxt.m_fibs.put( tag, 0 ); break;
        case 1 : ctxt.m_fibs.put( tag, 1 ); break;
        default : 
            // get previous 2 results
            tree_type f_1; ctxt.m_fibs.get( tag - 1, f_1 );
            tree_type f_2; ctxt.m_fibs.get( tag - 2, f_2 );
            // put our result
            ctxt.m_fibs.put( tag, f_1 + f_2 );
    }
    return CnC::CNC_Success;
}

int main( int argc, char* argv[] )
{
    int n = 100;
    // eval command line args
    if( argc < 2 ) {
        std::cerr << "usage: " << argv[0] << " n\nUsing default height " << n << std::endl;
    } else n = atol( argv[1] );

    // create context
    tree_context ctxt;

    // put tags to initiate evaluation
    for( int i = 0; i <= n; ++i ) ctxt.m_tags.put( i );

    // wait for completion
    ctxt.wait(); 

    // get result
    tree_type res2;
    ctxt.m_fibs.get( n, res2 );

    // print result
    std::cout << "tree (" << n << "): " << res2 << std::endl;

    return 0;
}
