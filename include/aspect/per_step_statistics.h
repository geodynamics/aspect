class PerStepStatistics
{
  public:
    void new_line ()
      {
	entries.push_back (std::list<std::pair<std::string,std::string> >());
      }

    void add_entry (const std::string &name,
		    const std::string &value)
      {
	entries.back()[name] = value;

	if (keys.find (name) == keys.end())
	  keys.push_back (name);
      }

    void add_entries (const std::map<std::string,std::string> &new_entries)
      {
	entries.back().insert (new_entries.begin(),
			       new_entries.end(),
			       entries.back().end());

	for (p = new_entries.begin(); p != new_entries.end(); ++p)
	  if (keys.find (p->first) == keys.end())
	    keys.push_back (p->first);
      }



                     /**
                      * Read or write the data of this object to or
                      * from a stream for the purpose of serialization
                      */
    template <class Archive>
    void serialize(Archive & ar, const unsigned int)
      {
	ar & entries & keys;
      }


  private:
    std::list<std::map<std::string,std::string> >
    entries;

    std::list<std::string> keys;
};
