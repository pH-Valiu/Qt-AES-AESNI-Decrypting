#include <QCoreApplication>
#include <QTest>
#include "gcm/gfmul.h"
#include "gcm/gcm.h"

#ifdef __cplusplus
#include "unit_test/aestest.h"
#endif

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    gfmul_test();
    //gcm_test();
    //AesTest test1;
    //return QTest::qExec(&test1);


    return 1;
}


