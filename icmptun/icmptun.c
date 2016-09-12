/*
 * Matteo Cicuttin (C) 2016 - IPoICMP
 *
 * IPoICMP kernel mode driver for FreeBSD...writing it just for fun.
 * Not complete and don't know when/if it will be finished.
 * Released under BSD License
 */

#include <sys/param.h>
#include <sys/proc.h>
#include <sys/module.h>
#include <sys/kernel.h>
#include <sys/systm.h>
#include <sys/socket.h>
#include <sys/sockio.h>
#include <sys/mbuf.h>
#include <sys/protosw.h>

#include <sys/conf.h>
#include <sys/uio.h>
#include <sys/malloc.h>

#include <sys/types.h>
#include <sys/sysctl.h>

#include <net/if.h>
#include <net/if_types.h>
#include <net/if_var.h>

#include <net/netisr.h>
#include <net/route.h>
#include <net/vnet.h>
#include <netinet/in.h>
#include <netinet/in_systm.h>
#include <netinet/ip.h>
#include <netinet/ip_icmp.h>
#include <netinet/ip_var.h>
#include <net/bpf.h>

static d_open_t     icmptun_open;
static d_close_t    icmptun_close;
static d_read_t     icmptun_read;
static d_write_t    icmptun_write;

static pr_input_t   icmp_input_hook;

extern struct protosw inetsw[];

static struct cdevsw icmptun_cdevsw = {
    .d_version =    D_VERSION,
    .d_open =       icmptun_open,
    .d_close =      icmptun_close,
    .d_read =       icmptun_read,
    .d_write =      icmptun_write,
    .d_name =       "icmptun"
};

struct icmptun_softc {
    struct cdev         *icmptun_dev;
    struct ifnet        *icmptun_ifp;
};


struct cdev         *icmptun_dev;

static int icmptundebug;
static int icmptun_rx_dump;
static int icmptun_tx_dump;

SYSCTL_INT(_debug, OID_AUTO, if_icmptun_debug, CTLFLAG_RW, &icmptundebug, 0, "");
SYSCTL_INT(_debug, OID_AUTO, if_icmptun_rx_dump, CTLFLAG_RW, &icmptun_rx_dump, 0, "");
SYSCTL_INT(_debug, OID_AUTO, if_icmptun_tx_dump, CTLFLAG_RW, &icmptun_tx_dump, 0, "");

static const char icmptunname[] = "icmptun";
static MALLOC_DEFINE(M_ICMPTUN, icmptunname, "IPoICMP Interface");

static struct icmptun_softc *isc;

static void
icmp_input_hook(struct mbuf *m, int off)
{
    struct icmp     *icp;
    int             hlen = off;
    
    m->m_len -= hlen;
    m->m_data += hlen;
    
    icp = mtod(m, struct icmp *);
    
    m->m_len += hlen;
    m->m_data -= hlen;
    
    m = m_pullup(m, m->m_len);
    
    if (icmptun_rx_dump)
    {
        printf("received icmp type %d code %d\n", icp->icmp_type, icp->icmp_code);
        hexdump(m->m_data, m->m_len, NULL, 0);
    }
    
    icmp_input(m, off);
}

static int
icmptun_open(struct cdev *dev, int oflags, int devtype, struct thread *td)
{
    return 0;
}

static int
icmptun_close(struct cdev *dev, int fflag, int devtype, struct thread *td)
{
    return 0;
}

static int
icmptun_read(struct cdev *dev, struct uio *uio, int ioflag)
{
    int error = 0;
    
    return error;
}

static int
icmptun_write(struct cdev *dev, struct uio *uio, int ioflag)
{
    int error = 0;
    
    return error;
}

static void
icmptun_ifstart(struct ifnet *ifp)
{
}

static void
icmptun_ifinit(void *xtp)
{
    struct icmptun_softc    *tp = (struct icmptun_softc *)xtp;
    struct ifnet            *ifp = tp->icmptun_ifp;
    
    icmptun_ifstart(ifp);
}

static int
icmptun_ifioctl(struct ifnet *ifp, u_long cmd, caddr_t data)
{
    struct ifreq *ifr = (struct ifreq *)data;
    int error = 0, mask;
    
    switch (cmd) {
        case SIOCSIFADDR:
            ifp->if_flags |= IFF_UP;
            ifp->if_drv_flags |= IFF_DRV_RUNNING;
            /*
             * Everything else is done at a higher level.
             */
            break;
            
        case SIOCADDMULTI:
        case SIOCDELMULTI:
            if (ifr == NULL) {
                error = EAFNOSUPPORT;		/* XXX */
                break;
            }
            switch (ifr->ifr_addr.sa_family) {
                    
#ifdef INET
                case AF_INET:
                    break;
#endif
#ifdef INET6
                case AF_INET6:
                    break;
#endif
                    
                default:
                    error = EAFNOSUPPORT;
                    break;
            }
            break;
            
        case SIOCSIFMTU:
            ifp->if_mtu = ifr->ifr_mtu;
            break;
            
        case SIOCSIFFLAGS:
            break;
            
        case SIOCSIFCAP:
            mask = ifp->if_capenable ^ ifr->ifr_reqcap;
            if ((mask & IFCAP_RXCSUM) != 0)
                ifp->if_capenable ^= IFCAP_RXCSUM;
            if ((mask & IFCAP_TXCSUM) != 0)
                ifp->if_capenable ^= IFCAP_TXCSUM;
            if ((mask & IFCAP_RXCSUM_IPV6) != 0) {
#if 0
                ifp->if_capenable ^= IFCAP_RXCSUM_IPV6;
#else
                error = EOPNOTSUPP;
                break;
#endif
            }
            if ((mask & IFCAP_TXCSUM_IPV6) != 0) {
#if 0
                ifp->if_capenable ^= IFCAP_TXCSUM_IPV6;
#else
                error = EOPNOTSUPP;
                break;
#endif
            }
            ifp->if_hwassist = 0;
            break;
            
        default:
            error = EINVAL;
    }
    return (error);
}

static int
icmptun_output(struct ifnet *ifp, struct mbuf *m, const struct sockaddr *dst,
               struct route *ro)
{
    u_int32_t       af;
    
    af = dst->sa_family;
    
    m = m_pullup(m, m->m_len);
    if (m == NULL)
        return ENOBUFS;
    
    char *data = mtod(m, char *);
    
    hexdump(data, m->m_len, NULL, 0);
    
    m_freem(m);
    return EHOSTUNREACH;
}


static struct icmptun_softc *
icmptun_create(const char *name, struct cdev *dev)
{
    struct icmptun_softc *sc;
    struct ifnet *ifp;
    
    sc = malloc(sizeof(*sc), M_ICMPTUN, M_WAITOK | M_ZERO);
    sc->icmptun_dev = dev;
    ifp = sc->icmptun_ifp = if_alloc(IFT_TUNNEL);
    if (ifp == NULL)
        panic("%s%d: if_alloc() failed\n", name, dev2unit(dev));
    
    if_initname(ifp, name, dev2unit(dev));
    //ifp->if_init = icmptun_ifinit;
    //ifp->if_start = icmptun_ifstart;
    ifp->if_ioctl = icmptun_ifioctl;
    ifp->if_output = icmptun_output;
    ifp->if_mtu = 1024;
    
    if_attach(ifp);
    
    bpfattach(ifp, DLT_NULL, sizeof(u_int32_t));
    
    return sc;
}

static void
icmptun_destroy(struct icmptun_softc *sc)
{
    struct cdev *dev;
    dev = sc->icmptun_dev;
    bpfdetach(sc->icmptun_ifp);
    if_detach(sc->icmptun_ifp);
    if_free(sc->icmptun_ifp);
    free(sc, M_ICMPTUN);
}

static int
icmptun_modevent(module_t mod __unused, int event, void *arg __unused)
{
    int error = 0;
    
    switch (event)
    {
        case MOD_LOAD:
            icmptun_dev = make_dev(&icmptun_cdevsw, 0, UID_ROOT, GID_WHEEL,
                                   0600, "icmptun");
            isc = icmptun_create(icmptunname, icmptun_dev);
            inetsw[ip_protox[IPPROTO_ICMP]].pr_input = icmp_input_hook;
            printf("ICMPTUN driver loaded\n");
            break;
    
        case MOD_UNLOAD:
            inetsw[ip_protox[IPPROTO_ICMP]].pr_input = icmp_input;
            icmptun_destroy(isc);
            destroy_dev(icmptun_dev);
            printf("ICMPTUN driver unloaded\n");
            break;
    
        default:
            error = EOPNOTSUPP;
            break;
    }
    
    return (error);
}

static moduledata_t icmptun_mod = {
    "icmptun",
    icmptun_modevent,
    NULL
};

DECLARE_MODULE(icmptun, icmptun_mod, SI_SUB_DRIVERS, SI_ORDER_MIDDLE);


