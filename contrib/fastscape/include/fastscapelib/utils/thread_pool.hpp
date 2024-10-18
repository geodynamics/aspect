#ifndef FASTSCAPELIB_UTILS_THREAD_POOL_HPP
#define FASTSCAPELIB_UTILS_THREAD_POOL_HPP

#include <atomic>
#include <chrono>
#include <functional>
#include <future>
#include <iostream>
#include <thread>
#include <condition_variable>


namespace fastscapelib
{
    template <class T>
    class thread_pool
    {
    public:
        using job_type = std::function<void()>;

        explicit thread_pool(size_t size);

        ~thread_pool();

        void set_tasks(std::vector<job_type>& jobs_);

        void run_tasks();

        void pause();

        void resume();

        bool paused() const;

        bool was_empty() const;

        void wait() const;

        void stop();

        std::size_t size() const;

        bool stopped() const;

        void start();

        bool started() const;

        void resize(std::size_t size);

        class [[nodiscard]] blocks
        {
        public:
            /**
             * @brief Construct a `blocks` object with the given specifications.
             *
             * @param first_index_ The first index in the range.
             * @param index_after_last_ The index after the last index in the range.
             * @param num_blocks_ The desired number of blocks to divide the range into.
             * @param min_size_ The minimum size of each blocks (ignored if zero, default).
             */
            blocks(const T& first_index_,
                   const T& index_after_last_,
                   const std::size_t num_blocks_,
                   const std::size_t min_size_ = 0);

            /**
             * @brief Get the first index of a block.
             *
             * @param block The block number.
             * @return The first index.
             */
            [[nodiscard]] T start(const std::size_t block) const;

            /**
             * @brief Get the index after the last index of a block.
             *
             * @param block The block number.
             * @return The index after the last index.
             */
            [[nodiscard]] T end(const std::size_t block) const;

            /**
             * @brief Get the number of blocks. Note that this may be different than the desired
             * number of blocks that was passed to the constructor.
             *
             * @return The number of blocks.
             */
            [[nodiscard]] std::size_t num_blocks() const;

        private:
            std::size_t m_block_size = 0;  ///< size of each block (except possibly the last block)
            T m_first_index = 0;           ///< first index in the range
            T m_index_after_last = 0;      ///< index after the last index in the range
            std::size_t m_num_blocks = 0;  /// < number of blocks
            std::size_t m_remainder
                = 0;  ///< remainder obtained after dividing the total size by the number of blocks

        };  // class blocks

        template <typename F>
        void run_blocks(const T first_index,
                        const T index_after_last,
                        F&& func,
                        const std::size_t min_size = 0);

    private:
        std::vector<std::thread> m_workers;
        std::vector<job_type>* p_jobs;
        std::vector<job_type> m_pause_jobs;
        std::atomic_bool m_stopped;
        std::vector<std::atomic<std::uint8_t>> m_has_job;
        std::size_t m_size;
        bool m_started = false, m_paused = false;
        std::atomic<std::size_t> m_paused_count = 0;

        std::condition_variable m_cv;
        std::mutex m_cv_m;

        void init_pause_jobs();
    };
}  // namespace fastscapelib

#include "./impl/thread_pool_inl.hpp"

#endif  // FASTSCAPELIB_UTILS_THREAD_POOL_HPP
