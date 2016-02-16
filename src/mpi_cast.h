#pragma once

#if defined(USE_MPI)

#include <string>
#include <cstring>
#include <vector>
#include <mpi.h>


namespace cmpi
{

  template<std::size_t N>
  void safe_string_to_charr (std::string const &str, char (&charr)[N])
  {
    std::size_t i(0U);
    for (; i<(N-1) && i<str.length(); ++i)
    {
      charr[i] = str[i];
    }
    charr[++i] = '\0';
  }

  void perform_checked (int err_code, std::string const & prefix="");

  class CMPI
  {
  protected:
    // mpi client and status
    MPI_Comm    m_client;
    MPI_Status  m_status;
    // connected flag?
    bool m_connected;
    // initialization checkvar
    static bool init_done;
    // instance counter for finalization
    static std::size_t instances;
  public:
    // construct empty mpi client
    CMPI () 
      : m_client(0), m_status(), m_connected(false) 
    { 
      ++CMPI::instances; 
      if (!CMPI::init_done) CMPI::init();
    }
    // construct new mpi client for given server
    CMPI (std::string server_name) 
      : m_client(0), m_status(), m_connected(false) 
    {
      ++CMPI::instances;
      if (!CMPI::init_done) CMPI::init();
      connect(server_name);
    }
    // connect to server
    void connect (std::string server_name);
    // destroy and finalize mpi if necessary
    ~CMPI ();
    // Status getter
    MPI_Status const & status () const { return m_status; }
    bool connected () const { return m_connected; }
    // mpi initialization
    static void init ();
    // send and recieve wrapper
    void send (void * data, int count, MPI_Datatype type, int dest, int tag) const;
    void recv (void * data, int count, MPI_Datatype type, int src, int tag);
    // send and recieve vectors, specialized in .cc
    template<class T>
    void send (std::vector<T> &data_vector, int dest, int tag) const;
    template<class T>
    void recv (std::vector<T> &data_vector, int dest, int tag);
    // send and recieve single values, specialized in .cc
    template<class T>
    void send (T & data, int dest, int tag) const;
    template<class T>
    void recv (T & data, int dest, int tag);


  };

  class CMPI_Srv : public CMPI
  {
  private:
    CMPI_Srv & operator= (CMPI_Srv const &);
  protected:
    char m_port_name[MPI_MAX_PORT_NAME];
    std::string m_publish_name;
    bool m_open, m_published, m_accepting;

    void open();
    void publish();
    void accept();
  public:
    CMPI_Srv()
      : CMPI(), m_port_name(), m_publish_name(),
      m_open(false), m_published(false), m_accepting(false)
    { }
    CMPI_Srv(std::string const & server_publish_name)
      : CMPI(), m_port_name(), m_publish_name(server_publish_name),
      m_open(false), m_published(false), m_accepting(false)
    { open_publish_accept(server_publish_name); }
    ~CMPI_Srv();
    void open_publish_accept (std::string const & publish_name)
    {
      m_publish_name = publish_name;
      open();
      publish();
      accept();
    }
  };


}

#endif

