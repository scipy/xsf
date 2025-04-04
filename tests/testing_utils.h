#pragma once

#include <algorithm>
#include <complex>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>

#include <arrow/io/file.h>
#include <parquet/stream_reader.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/generators/catch_generators_adapters.hpp>

#include <xsf/fp_error_metrics.h>

namespace {

template <typename T>
struct remove_complex {
    using type = T;
};

template <typename T>
struct remove_complex<std::complex<T>> {
    using type = T;
};

template <typename T>
using remove_complex_t = typename remove_complex<T>::type;

template <typename T>
class TableReader {
  public:
    TableReader(const std::string &file_path) {
        PARQUET_ASSIGN_OR_THROW(infile_, arrow::io::ReadableFile::Open(file_path));
        parquet::StreamReader stream{parquet::ParquetFileReader::Open(infile_)};
        stream_ = std::make_unique<parquet::StreamReader>(std::move(stream));
        file_path_ = std::move(file_path);
    }

    T next() {
        T row;
        fill_row(row);
        return row;
    }

    bool eof() const { return stream_->eof(); }

  private:
    void fill_row(T &elements) {
        if constexpr (std::is_scalar_v<remove_complex_t<T>>) {
            fill_element(elements);
        } else {
            std::apply([this](auto &...x) { (fill_element(x), ...); }, elements);
        }
        stream_->EndRow();
    }

    template <typename U>
    void fill_element(U &element) {
        if constexpr (std::is_same_v<U, std::complex<remove_complex_t<U>>>) {
            using V = remove_complex_t<U>;
            V real;
            V imag;
            *stream_ >> real >> imag;
            element = U(real, imag);
        } else {
            *stream_ >> element;
        }
    }

    std::shared_ptr<arrow::io::ReadableFile> infile_;
    std::unique_ptr<parquet::StreamReader> stream_;
    std::string file_path_;
};

template <typename T1, typename T2, typename T3>
class XsfTestCaseGenerator final : public Catch::Generators::IGenerator<std::tuple<T1, T2, T3>> {
  public:
    XsfTestCaseGenerator(
        std::unique_ptr<TableReader<T1>> input_reader, std::unique_ptr<TableReader<T2>> output_reader,
        std::unique_ptr<TableReader<T3>> tol_reader
    )
        : input_reader_(std::move(input_reader)), output_reader_(std::move(output_reader)),
          tol_reader_(std::move(tol_reader)) {
        if (!next()) {
            throw std::runtime_error("XsfTestCaseGenerator received an empty table\n");
        }
    }

    std::tuple<T1, T2, T3> const &get() const override { return current_case_; }

    bool next() override {
        if (input_reader_->eof() || output_reader_->eof() || tol_reader_->eof()) {
            return false;
        }
        current_case_ = std::make_tuple(input_reader_->next(), output_reader_->next(), tol_reader_->next());
        return true;
    }

  private:
    std::unique_ptr<TableReader<T1>> input_reader_;
    std::unique_ptr<TableReader<T2>> output_reader_;
    std::unique_ptr<TableReader<T3>> tol_reader_;
    std::tuple<T1, T2, T3> current_case_;
};

template <typename T1, typename T2, typename T3>
Catch::Generators::GeneratorWrapper<std::tuple<T1, T2, T3>> xsf_test_cases(
    const std::filesystem::path &input_path, const std::filesystem::path &output_path,
    const std::filesystem::path &tol_path
) {
    auto input_reader = std::make_unique<TableReader<T1>>(input_path.string());
    auto output_reader = std::make_unique<TableReader<T2>>(output_path.string());
    auto tol_reader = std::make_unique<TableReader<T3>>(tol_path.string());

    return Catch::Generators::GeneratorWrapper<std::tuple<T1, T2, T3>>(
        Catch::Detail::make_unique<XsfTestCaseGenerator<T1, T2, T3>>(
            std::move(input_reader), std::move(output_reader), std::move(tol_reader)
        )
    );
}


template <typename T>
T adjust_tolerance(T tol) {
    // Add some wiggle room to tolerance from table.
    return 2 * std::max(std::numeric_limits<T>::epsilon(), tol);
}


std::string get_platform_str() {
    /* This is going to get a string "<compiler>-<os>-<architecture>" using conditional
     * compilation, but for now we're just stubbing things out. */
    return "gcc-linux-x86_64";
}



} // namespace


#define SET_FP_FORMAT()							\
    Catch::StringMaker<double>::precision = std::numeric_limits<double>::max_digits10; \
    Catch::StringMaker<float>::precision = std::numeric_limits<float>::max_digits10;
