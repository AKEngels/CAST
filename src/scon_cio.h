#if !defined(SCON_IO_HEADER)
#define SCON_IO_HEADER

#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

#pragma warning(disable : 4996)

namespace scon
{

  namespace io
  {
    template<class T>
    std::vector<char> char_buffer_from_vector(std::vector<T> const & vec)
    {
      std::size_t const
        nv = vec.size(),
        nb = nv*sizeof(typename std::vector<T>::value_type),
        nc = nb / sizeof(char);
      std::vector<char> buf(nc, 0);
      std::memcpy(buf.data(), vec.data(), nb);
      return buf;
    }

    class cfile
    {
    protected:
      FILE * m_file;
      bool open(char const * const filename, char const * const mode)
      {
        cfile tmp(filename, mode);
        std::swap(m_file, tmp.m_file);
        return m_file != NULL;
      }
      cfile() : m_file(NULL) { }
      cfile(char const * const filename, char const * const mode)
        : m_file(fopen(filename, mode))
      { }

      bool is_open() const { return m_file != NULL; }
    public:
      bool close()
      {
        if (is_open())
        {
          int const ret = fclose(m_file);
          m_file = NULL;
          return ret == 0;
        }
        return true;
      }
      ~cfile()
      {
        close();
      }
    };

    template<bool out = false>
    class binfile
      : public cfile
    {
    public:
      binfile() : cfile() { }
      binfile(char const * const filename)
        : cfile(filename, "rb")
      { }
      bool open(char const * const filename)
      {
        binfile tmp(filename);
        std::swap(m_file, tmp.m_file);
        return m_file != NULL;
      }
      bool read(void * const ptr, std::size_t const element_size, std::size_t const element_count)
      {
        if (is_open() && !out)
        {
          std::size_t const num_read = fread(ptr, element_size, element_count, m_file);
          if (num_read != element_count)
          {
            throw std::system_error(std::make_error_code(std::io_errc::stream), "Read from file failed.");
          }
          return true;
        }
        return false;
      }

      bool read(void * const ptr, std::size_t const bytes)
      {
        if (is_open() && !out)
        {
          std::size_t const num_read = fread(ptr, bytes, 1, m_file);
          if (num_read == 0)
          {
            throw std::system_error(std::make_error_code(std::io_errc::stream), "Read from file failed.");
          }
          return true;
        }
        return false;
      }
      std::vector<char> read()
      {
        if (is_open() && !out)
        {
          fseek(m_file, 0, SEEK_END);
          long filesize = ftell(m_file);
          fseek(m_file, 0, SEEK_SET);
          if (filesize > 0)
          {
            std::vector<char> content(static_cast<std::size_t(filesize));
            read(content.data(), filesize);
            return content;
          }
          else if (filesize == -1L)
          {
            throw std::system_error(std::make_error_code(std::io_errc::stream), "Read from file failed.");
          }
        }
        return std::vector<char>();
      }
    };

    template<>
    class binfile < true > 
      : public cfile
    {
    public:
      binfile() : cfile() { }
      binfile(char const * const filename)
        : cfile(filename, "wb")
      { }
      bool open(char const * const filename)
      {
        binfile tmp(filename);
        std::swap(m_file, tmp.m_file);
        return m_file != NULL;
      }
      bool write(void const * const ptr, std::size_t const element_size, std::size_t const element_count)
      {
        if (is_open())
        {
          std::size_t const written = fwrite(ptr, element_size, element_count, m_file);
          if (written != element_count)
          {
            throw std::system_error(std::make_error_code(std::io_errc::stream), "Write to file failed.");
          }
          return true;
        }
        return false;
      }

      bool write(void const * const ptr, std::size_t const bytes)
      {
        if (is_open())
        {
          std::size_t const written = fwrite(ptr, bytes, 1, m_file);
          if (written == 0)
          {
            throw std::system_error(std::make_error_code(std::io_errc::stream), "Write to file failed.");
          }
          return true;
        }
        return false;
      }

      bool write(std::vector<char> const &buffer)
      {
        return write(buffer.data(), buffer.size()*sizeof(char));
      }

    };

    typedef binfile<false> bin_in_file;
    typedef binfile<true> bin_out_file;

  }
  
}

#endif
