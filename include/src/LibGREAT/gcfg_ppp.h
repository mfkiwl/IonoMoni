#ifndef GCFG_PPP_H
#define GCFG_PPP_H

#include <string>
#include <iostream>
#include <signal.h>

#include "gall/gallprec.h"
#include "gall/gallpcv.h"
#include "gall/gallobs.h"
#include "gall/gallotl.h"
#include "gall/gallbias.h"
#include "gproc/gpppflt.h"
#include "gproc/gpreproc.h"
#include "gio/gio.h"
#include "gio/gfile.h"
#include "gio/grtlog.h"       // ✅ 新增：日志模块
#include "spdlog/spdlog.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"
#include "gdata/gdata.h"
#include "gcoders/gcoder.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexc.h"
#include "gcoders/rinexn.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/sp3.h"
#include "gcoders/atx.h"
#include "gcoders/blq.h"
#include "gmodels/gpcv.h"
#include "gmodels/gbancroft.h"
#include "gset/gsetgen.h"
#include "gset/gsetinp.h"
#include "gset/gsetout.h"
#include "gset/gsetproc.h"
#include "gset/gsetgnss.h"
#include "gset/gsetflt.h"
#include "gset/gsetrec.h"
#include "gproc/gpvtflt.h"
#include "gdata/gifcb.h"
#include "gcoders/upd.h"
#include "gcoders/sp3.h"
#include "gcoders/rinexo.h"
#include "gcoders/rinexn.h"
#include "gcoders/rinexc.h"
#include "gcoders/poleut1.h"
#include "gcoders/ifcb.h"
#include "gcoders/gcoder.h"
#include "gcoders/dvpteph405.h"
#include "gcoders/blq.h"
#include "gcoders/biasinex.h"
#include "gcoders/biabernese.h"
#include "gcoders/atx.h"
#include "gdata/gnavde.h"

using namespace std;
using namespace pugi;

namespace gnut
{

    class t_gcfg_ppp : public t_gsetgen,
        public t_gsetinp,
        public t_gsetout,
        public t_gsetgnss,
        public t_gsetproc,
        public t_gsetflt,
        public t_gsetrec,
        public t_gsetamb
    {
    public:
        /** @brief default constructor. */
        t_gcfg_ppp();

        /** @brief default destructor. */
        ~t_gcfg_ppp();

        /** @brief settings check. */
        void check();

        /** @brief settings help. */
        void help();

        /** @brief 读取XML自定义字段，比如<function>XXX</function> */
        string comment(const string& key)
        {
            _gmutex.lock();
            string val = _doc.child("config").child("gen").child(key.c_str()).child_value();
            _gmutex.unlock();
            return trim(val);
        }

        /** @brief 支持直接给XML路径（不是argc/argv方式） */
        void arg(const string& xmlfile, bool xml = true, bool log = true);
        void arg(int argc, char** argv, bool xml, bool log);

        // 🔥【新增】日志控制函数
        void init_logger();  // 初始化内部日志器
        std::shared_ptr<spdlog::logger> logger();  // 提供日志器接口给外部用

    protected:
    private:
        t_grtlog _grtlog;   // 🔥【新增】Great标准日志对象
    };

} // namespace gnut

#endif
