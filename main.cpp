#include <QCoreApplication>
#include <QTest>
#include "gcm/gfmul.h"

#ifdef __cplusplus
#include "unit_test/aestest.h"
#endif

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    //AesTest test1;
    //return QTest::qExec(&test1);


    gfmul_test();
    return 1;
}


