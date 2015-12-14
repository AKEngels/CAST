#include "mpi_cast.h"

#if defined(USE_MPI)

#include <iostream>
#include <climits>
#include <stdexcept>


bool cmpi::CMPI::init_done = false;

std::size_t cmpi::CMPI::instances = 0U;

void cmpi::perform_checked (int err_code, std::string const & prefix)
{
  if (err_code != MPI_SUCCESS)
  {
    int error_len(0);
    char error_string[MPI_MAX_ERROR_STRING]; 
    MPI_Error_string(err_code, error_string, &error_len);
    throw std::runtime_error((prefix + std::string(error_string)).c_str());
  } 
  //else
  //{
  //  std::cout << prefix << " success." << std::endl;
  //}
}

void cmpi::CMPI::init (void) 
{
  int rank; 
  // Try to initialize MPI
  perform_checked(MPI_Init( 0, NULL ), "ERR_MPI_Init: ");
  // try to Set proper Comm Size
  perform_checked(MPI_Comm_size(MPI_COMM_WORLD, &rank), "ERR_MPI_Comm_size: ");
  std::cout << "Rank is: " << rank << std::endl;
  CMPI::init_done = true;
}

void cmpi::CMPI::connect (std::string server_name)
{
  char port_name[MPI_MAX_PORT_NAME], srv_name[256] = ""; 
  safe_string_to_charr(server_name, srv_name);
  std::cout << "MPI::Server is: " << server_name << " which is " << srv_name << std::endl;
  // Lookup Server port
  perform_checked(MPI_Lookup_name(srv_name, MPI_INFO_NULL, port_name), "ERR_MPI_Lookup_name : ");
  std::cout << "MPI::Server found (" << port_name << ") Connecting... " << std::endl;
  // Try to connect
  perform_checked(MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &m_client), "ERR_MPI_Comm_connect : ");
  m_connected = true;
  std::cout << "MPI::Connected. " << std::endl;
}

cmpi::CMPI::~CMPI (void)
{
  if (connected()) MPI_Comm_disconnect(&m_client);
  --CMPI::instances;
  if (CMPI::init_done && CMPI::instances == 0) MPI_Finalize();
}

void cmpi::CMPI::send (void const * data, int size, MPI_Datatype type, int dest, int tag) const
{
  perform_checked(MPI_Send(data, size, type, dest, tag, m_client), "ERR_MPI_Send : ");
}

void cmpi::CMPI::recv (void * data, int size, MPI_Datatype type, int src, int tag)
{
  perform_checked(MPI_Recv(data, size, type, src, tag, m_client, &m_status), "ERR_MPI_Recv : ");
}

namespace cmpi
{

  /*
  SEND AND RECIEVE VECTORS
  */


  //! double

  template<>
  void CMPI::send<double> (std::vector<double> &data_vector, int dest, int tag) const
  {
    if (data_vector.size() > INT_MAX) throw std::runtime_error("ERR_MPI_SEND_DOUBLE_VECTOR_TOO_LARGE");
    send(data_vector.data(), static_cast<int>(data_vector.size()), MPI_DOUBLE, dest, tag);
  }
  template<>
  void CMPI::recv<double> (std::vector<double> &data_vector, int src, int tag)
  {
    if (data_vector.size() > INT_MAX) throw std::runtime_error("ERR_MPI_RECV_DOUBLE_VECTOR_TOO_LARGE");
    recv(data_vector.data(), static_cast<int>(data_vector.size()), MPI_DOUBLE, src, tag);
  }

  //! float

  template<>
  void CMPI::send<float> (std::vector<float> &data_vector, int dest, int tag) const
  {
    if (data_vector.size() > INT_MAX) throw std::runtime_error("ERR_MPI_SEND_FLOAT_VECTOR_TOO_LARGE");
    send(data_vector.data(), static_cast<int>(data_vector.size()), MPI_FLOAT, dest, tag);
  }
  template<>
  void cmpi::CMPI::recv<float> (std::vector<float> &data_vector, int src, int tag)
  {
    if (data_vector.size() > INT_MAX) throw std::runtime_error("ERR_MPI_RECV_FLOAT_VECTOR_TOO_LARGE");
    recv(data_vector.data(), static_cast<int>(data_vector.size()), MPI_FLOAT, src, tag);
  }

  //! int

  template<>
  void CMPI::send<int> (std::vector<int> &data_vector, int dest, int tag) const
  {
    if (data_vector.size() > long long(INT_MAX)) throw std::runtime_error("ERR_MPI_SEND_INT_VECTOR_TOO_LARGE");
    send(data_vector.data(), static_cast<int>(data_vector.size()), MPI_INT, dest, tag);
  }
  template<>
  void CMPI::recv<int> (std::vector<int> &data_vector, int src, int tag)
  {
    if (data_vector.size() > INT_MAX) throw std::runtime_error("ERR_MPI_RECV_INT_VECTOR_TOO_LARGE");
    recv(data_vector.data(), static_cast<int>(data_vector.size()), MPI_INT, src, tag);
  }

  //! char

  template<>
  void cmpi::CMPI::send<char> (std::vector<char> &data_vector, int dest, int tag) const
  {
    if (data_vector.size() > INT_MAX) throw std::runtime_error("ERR_MPI_SEND_CHAR_VECTOR_TOO_LARGE");
    send(data_vector.data(), static_cast<int>(data_vector.size()), MPI_CHAR, dest, tag);
  }
  template<>
  void cmpi::CMPI::recv<char> (std::vector<char> &data_vector, int src, int tag)
  {
    if (data_vector.size() > INT_MAX) throw std::runtime_error("ERR_MPI_RECV_CHAR_VECTOR_TOO_LARGE");
    recv(data_vector.data(), static_cast<int>(data_vector.size()), MPI_CHAR, src, tag);
  }

  /*
  SEND AND RECIEVE SINGLE VALUES
  */

  //! double
  template<> void cmpi::CMPI::send<double> (double & data, int dest, int tag) const { send(&data, 1, MPI_DOUBLE, dest, tag); }
  template<> void cmpi::CMPI::recv<double> (double & data, int src, int tag) { recv(&data, 1, MPI_DOUBLE, src, tag); }

  //! float
  template<> void cmpi::CMPI::send<float> (float & data, int dest, int tag) const { send(&data, 1, MPI_FLOAT, dest, tag); }
  template<> void cmpi::CMPI::recv<float> (float & data, int src, int tag) { recv(&data, 1, MPI_FLOAT, src, tag); }

  //! int
  template<> void cmpi::CMPI::send<int> (int & data, int dest, int tag) const { send(&data, 1, MPI_INT, dest, tag); }
  template<> void cmpi::CMPI::recv<int> (int & data, int src, int tag) { recv(&data, 1, MPI_INT, src, tag); }

  //! char
  template<> void cmpi::CMPI::send<char> (char & data, int dest, int tag) const { send(&data, 1, MPI_CHAR, dest, tag); }
  template<> void cmpi::CMPI::recv<char> (char & data, int src, int tag) { recv(&data, 1, MPI_CHAR, src, tag); }

  //! string
  template<> void cmpi::CMPI::send<std::string> (std::string & data, int const dest, int const tag) const 
  {
    if (data.size() > INT_MAX) throw std::runtime_error("ERR_MPI_SEND_STRING_TOO_LARGE");
    send(static_cast<void*>(&data[0]), static_cast<int>(data.size()), MPI_CHAR, dest, tag);
  }
  template<> void cmpi::CMPI::recv<std::string>(std::string & data, int const src, int const tag)
  {
    if (data.size() > INT_MAX) throw std::runtime_error("ERR_MPI_RECV_STRING_TOO_LARGE");
    recv(static_cast<void*>(&data[0]), static_cast<int>(data.size()), MPI_CHAR, src, tag);
  }

}

cmpi::CMPI_Srv::~CMPI_Srv()
{
  if (m_accepting) 
    perform_checked(MPI_Comm_free(&m_client), "ERR_MPI_Comm_free : ");
  if (m_published) 
    perform_checked(MPI_Unpublish_name(&m_publish_name[0], MPI_INFO_NULL, m_port_name), "ERR_MPI_Unpublish_name : ");
  if (m_open) 
    perform_checked(MPI_Close_port(m_port_name), "ERR_MPI_Close_port : ");
}

void cmpi::CMPI_Srv::open()
{
  if (!m_open)
  {
    if (!CMPI::init_done) init();
    // Open MPI Port
    perform_checked(MPI_Open_port(MPI_INFO_NULL, m_port_name), "ERR_MPI_Open_port : ");
    m_open = true;
  }
}

void cmpi::CMPI_Srv::publish()
{
  if (!m_published && m_open)
  {
    m_publish_name.push_back('\0');
    // Open MPI Port
    perform_checked(MPI_Publish_name(&m_publish_name[0], MPI_INFO_NULL, m_port_name), "ERR_MPI_Publish_name : ");
    m_published = true;
  }
}

void cmpi::CMPI_Srv::accept()
{
  if (!m_accepting && m_open && m_published)
  {
    // Open MPI Port
    perform_checked(MPI_Comm_accept(m_port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &m_client), "ERR_MPI_Comm_accept : ");
    m_accepting = true;
  }    
}

#endif
