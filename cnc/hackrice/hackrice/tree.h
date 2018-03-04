#ifndef fib_H_ALREADY_INCLUDED
#define fib_H_ALREADY_INCLUDED

#include <cnc/cnc.h>
#include <cnc/debug.h>

// Forward declaration of the context class (also known as graph)
struct tree_context;

// The step classes

struct tree_step
{
	int execute(const int & t, tree_context & c) const;
};

// The context class
struct tree_context : public CnC::context< tree_context >
{
	// step collections
	CnC::step_collection< tree_step > m_steps;
	// Item collections
	CnC::item_collection< int, fib_type > m_tree;
	// Tag collections
	CnC::tag_collection< int > m_tags;

	// The context class constructor
	fib_context()
		: CnC::context< fib_context >(),
		// Initialize each step collection
		m_steps(*this),
		// Initialize each item collection
		m_trees(*this),
		// Initialize each tag collection
		m_tags(*this)
	{
		// Prescriptive relations
		m_tags.prescribes(m_steps, *this);
		// Consumer relations
		m_steps.consumes(m_fibs);
		// Producer relations
		m_steps.produces(m_fibs);
	}
};

#endif // fib_H_ALREADY_INCLUDED
