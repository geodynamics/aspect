/**
 * @file
 * @brief Provides implementation for the classical Union-Find data structure
 *
 * @author Guillaume Cordonnier
 */

#pragma once

#include "utils.hpp"
#include <vector>


namespace fastscapelib
{

    namespace detail
    {

        /**
         * @class union_find
         * @brief Union-Find data structure
         *
         * The union find is a standard data structure for storing and merging equivalence
         * classes. Initialy (after construction, or call to clear), the Union Find stores
         * size() distincts classes.
         * The classes can be merged with a call to merge(), and obtained through the call
         * to find().
         * Amortized complexity of Union and Find: optimal O(alpha(m, n))
         */
        template <class T>
        class union_find
        {
        public:
            /**
             * @brief UnionFind Constructor
             */
            union_find()
            {
            }

            /**
             * @brief UnionFind Constructor
             * @param _size the number of classes
             * Complexity O(_size)
             */
            union_find(size_t _size)
            {
                resize(_size);
            }

            /**
             * @brief clear
             * Restore the union-find data structure. The structure holds m_size different size()
             * distinct classes
             * Complexity O(size())
             */
            void clear()
            {
                size_t old_size = size();
                parent.clear();
                rank.clear();
                resize(old_size);
            }

            /**
             * @brief reserve
             * Allocate some memory for a given number of class
             * @param _size the number of classes
             */
            void reserve(size_t _size)
            {
                parent.reserve(_size);
                rank.reserve(_size);
            }

            /**
             * @brief push_back
             * append a new item at the end of the union find structure
             * @param c the class of the new item
             */

            void push_back(T c)
            {
                parent.push_back(c);
                rank.push_back(0);
                if (c != parent.size() - 1 && !rank[c])
                    rank[c] = 1;
            }

            /**
             * @brief resize
             * Resize the internal containers of union find. This only add classes if the new
             * size is larger. Else, some elements are removed, but it is possible that the
             * class returned by find is bigger than the new size.
             * Complexity O(_size - size())
             * @param _size the new size
             */
            void resize(size_t _size)
            {
                // size_t old_size = size();

                parent.resize(_size);
                rank.resize(_size, 0);

                // TODO: this causes seg fault when used from the py bindings
                // -> regenerate the whole sequence as a workaround
                // std::iota(parent.begin() + old_size, parent.end(), old_size);
                std::iota(parent.begin(), parent.end(), 0);
            }

            /**
             * @brief size
             * @return the initial number of elements in the union find structure
             */
            size_t size()
            {
                return parent.size();
            }

            /**
             * @brief merge two equivalence class
             * The new equivelence class can be represented either by x or y
             * @param x class to be merged
             * @param y class to be merged
             */
            void merge(T x, T y)
            {
                x = find(x);
                y = find(y);

                if (x != y)
                {
                    if (rank[x] < rank[y])
                        parent[x] = y;
                    else
                    {
                        parent[y] = x;
                        if (rank[x] == rank[y])
                            rank[x] += 1;
                    }
                }
            }

            /**
             * @brief find the equivalence class of an element
             * @param x the element for which the class is needed
             * @return the class representative
             */
            T find(T x)
            {
                // find class
                T c = x;
                while (c != parent[c])
                    c = parent[c];

                // set class
                while (x != parent[x])
                {
                    T t = parent[x];
                    parent[x] = c;
                    x = t;
                }

                // return class
                return c;
            }

        private:
            std::vector<T> parent;
            std::vector<size_t> rank;
        };
    }
}
