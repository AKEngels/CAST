#ifndef CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_MATRIX_ABSTRACTMATRIX_H_
#define CAST_INTERNALCOORDINATES_INTERNALCOORDINATES_MATRIX_ABSTRACTMATRIX_H_

#include<memory>

#include"../InternalCoordinatesAliases.h"

namespace internals{
    class AbstractMatrix {
    public:
        virtual float_type operator()(std::size_t const, std::size_t const) = 0;
        virtual AbstractMatrix & add(AbstractMatrix const&) = 0;
        virtual AbstractMatrix & subtract(AbstractMatrix const&) = 0;
        virtual AbstractMatrix & multiply(AbstractMatrix const&) = 0;
        virtual std::unique_ptr<AbstractMatrix> copy() const = 0;
        virtual std::size_t const cols() const = 0;
        virtual std::size_t const rows() const = 0;

    };
}

#endif