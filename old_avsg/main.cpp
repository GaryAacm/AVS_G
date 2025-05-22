#include "CommandParse.h"
#include "Context.h"

int main(int argc, char **argv)
{   
    time_t start_t, end_t;
    start_t = time(0);
    char * st = ctime(&start_t);
    cout << "start time: " << st << endl;

    Param *paramPtr = Param::GetInstance();
    CommandParse cmdParse;
    cmdParse.parseOptFromCmd(argc, argv);
    Context *contextPtr = new Context;

    switch (paramPtr->m_actionType)
    {
    case ACTIONTYPE_DOINDEX:
        contextPtr->doBuildIndex();
        break;
    case ACTIONTYPE_DOENCODE:
        contextPtr->initEncodeContext();
        break;
    case ACTIONTYPE_DODECODE:
        contextPtr->initDecodeContext();
        break;
    case ACTIONTYPE_EXTRACT:
        contextPtr->initExtractContext();
    default:
        break;
    }

    
    end_t = time(0);
    char * et = ctime(&end_t);
    cout << "end   time: " << et << endl;
    cout << "cost  time: " << (double)(end_t - start_t) << "s" << endl;

    RELEASEPTR(contextPtr);
    return 0;
}