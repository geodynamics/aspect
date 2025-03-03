#ifndef FASTSCAPELIB_UTILS_IMPL_THREAD_POOL_INL_HPP
#define FASTSCAPELIB_UTILS_IMPL_THREAD_POOL_INL_HPP


namespace fastscapelib
{
    template <class T>
    thread_pool<T>::thread_pool(std::size_t size)
        : m_stopped(false)
        , m_has_job(size)
        , m_size(size)
    {
        m_workers.reserve(size);
        init_pause_jobs();
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    thread_pool<T>::~thread_pool()
    {
        stop();
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    void thread_pool<T>::init_pause_jobs()
    {
        m_pause_jobs.resize(m_size);
        for (std::size_t i = 0; i < m_size; ++i)
            m_pause_jobs[i] = [this, i]()
            {
                std::unique_lock<std::mutex> lk(m_cv_m);
                ++m_paused_count;
                m_cv.wait(lk);
                --m_paused_count;
            };
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    void thread_pool<T>::set_tasks(std::vector<job_type>& jobs_)
    {
        p_jobs = &jobs_;
        assert(m_size == jobs_.size());
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    void thread_pool<T>::run_tasks()
    {
        if (!m_started)
            start();

        if (m_paused)
            resume();

        for (std::size_t i = 0; i < m_size; ++i)
            if ((*p_jobs)[i] != nullptr)
                m_has_job[i].store(1, std::memory_order_relaxed);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    void thread_pool<T>::pause()
    {
        if (!m_paused)
        {
            wait();
            set_tasks(m_pause_jobs);
            run_tasks();
            m_paused = true;

            while (m_paused_count != m_size)
            {
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    void thread_pool<T>::resume()
    {
        if (m_paused)
        {
            m_cv.notify_all();
            m_paused = false;
            wait();
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    bool thread_pool<T>::paused() const
    {
        return m_paused;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    bool thread_pool<T>::was_empty() const
    {
        for (std::size_t i = 0; i < m_size; ++i)
        {
            if (m_has_job[i].load(std::memory_order_relaxed))
                return false;
        }
        return true;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    void thread_pool<T>::wait() const
    {
        while (!was_empty())
        {
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    void thread_pool<T>::stop()
    {
        if (!m_stopped)
        {
            m_stopped = true;

            if (m_paused)
                resume();

            for (std::thread& worker : m_workers)
                worker.join();
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    std::size_t thread_pool<T>::size() const
    {
        return m_size;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    bool thread_pool<T>::stopped() const
    {
        return m_stopped;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    void thread_pool<T>::start()
    {
        if (!m_started)
        {
            m_started = true;

            for (size_t i = 0; i < m_size; ++i)
            {
                m_workers.emplace_back(
                    [this, i]
                    {
                        while (!m_stopped.load(std::memory_order_relaxed))
                        {
                            if (m_has_job[i].load(std::memory_order_relaxed))
                            {
                                (*p_jobs)[i]();
                                m_has_job[i].store(0, std::memory_order_relaxed);
                            }
                        }
                    });
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    bool thread_pool<T>::started() const
    {
        return m_started;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    void thread_pool<T>::resize(std::size_t size)
    {
        if (size != m_size)
        {
            m_size = size;
            stop();
            m_stopped = false;
            m_workers.clear();
            m_workers.reserve(size);
            m_has_job = std::vector<std::atomic<std::uint8_t>>(size);
            init_pause_jobs();
            m_started = false;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    template <typename F>
    void thread_pool<T>::run_blocks(const T first_index,
                                    const T index_after_last,
                                    F&& func,
                                    const std::size_t min_size)
    {
        std::vector<std::function<void()>> p_jobs(m_size);

        if (index_after_last > first_index)
        {
            const blocks blks(first_index, index_after_last, m_size, min_size);

            for (T i = 0; i < m_size; ++i)
            {
                if (i < blks.num_blocks())
                    p_jobs[i] = [i,
                                 func = std::forward<F>(func),
                                 start = blks.start(i),
                                 end = blks.end(i)]() { func(i, start, end); };
                else
                    p_jobs[i] = nullptr;
            }
            set_tasks(p_jobs);
            run_tasks();

            wait();
        };
    }

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    thread_pool<T>::blocks::blocks(const T& first_index_,
                                   const T& index_after_last_,
                                   const std::size_t num_blocks_,
                                   const std::size_t min_size_)
        : m_first_index(first_index_)
        , m_index_after_last(index_after_last_)
        , m_num_blocks(num_blocks_)
    {
        if (m_index_after_last > m_first_index)
        {
            const std::size_t total_size = static_cast<size_t>(m_index_after_last - m_first_index);
            if (m_num_blocks > total_size)
                m_num_blocks = total_size;

            if (total_size / m_num_blocks < min_size_)
                m_num_blocks = std::max(std::size_t{ 1 }, total_size / min_size_);

            m_block_size = total_size / m_num_blocks;
            m_remainder = total_size % m_num_blocks;
            if (m_block_size == 0)
            {
                m_block_size = 1;
                m_num_blocks = (total_size > 1) ? total_size : 1;
            }
        }
        else
        {
            m_num_blocks = 0;
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    T thread_pool<T>::blocks::start(const std::size_t block) const
    {
        return m_first_index + static_cast<T>(block * m_block_size)
               + static_cast<T>(block < m_remainder ? block : m_remainder);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    T thread_pool<T>::blocks::end(const std::size_t block) const
    {
        return (block == m_num_blocks - 1) ? m_index_after_last : start(block + 1);
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    template <class T>
    std::size_t thread_pool<T>::blocks::num_blocks() const
    {
        return m_num_blocks;
    }

}  // namespace fastscapelib

#endif  // FASTSCAPELIB_UTILS_IMPL_THREAD_POOL_INL_HPP
