#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        Vector2D interval = (end-start)/(num_nodes-1);
        for(int i=0; i<num_nodes; i++){
            Vector2D pos = Vector2D(start + i*interval);
            masses.push_back(new Mass(pos, node_mass, false));
        }
        for(int i=1; i<num_nodes; i++){
            springs.push_back(new Spring(masses[i-1], masses[i], k));
        }
//        Comment-in this part when you implement the constructor
       for (auto &i : pinned_nodes) {
           masses[i]->pinned = true;
       }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            auto dif = s->m2->position - s->m1->position;
            Vector2D f = s->k * dif.unit() * (dif.norm() - s->rest_length);
            s->m1->forces += f;
            s->m2->forces -= f;
        }

        
        for (auto &m : masses)
        {
            double damping_factor = 0.005;
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                // TODO (Part 2): Add global damping
                m->forces += m->mass * gravity;
                m->forces -= damping_factor * m->velocity;
                Vector2D a = m->forces / m->mass;               

                m->velocity += a * delta_t;
                m->position += m->velocity * delta_t;
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            auto a = s->m1->position;
            auto b = s->m2->position;
            auto l = (b - a).norm();
            auto half = (l - s->rest_length)/2;

            if(!s->m1->pinned) s->m1->position += half * (b-a).unit();
            if(!s->m2->pinned) s->m2->position -= half * (b-a).unit();
            
        }

        for (auto &m : masses)
        {   
            double damping_factor = 0.00005;
            if (!m->pinned)
            {
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass
                m->position = m->position + (1-damping_factor)*(m->position - m->last_position) + gravity*delta_t*delta_t;                             
                m->last_position = temp_position;
                // TODO (Part 4): Add global Verlet damping
            }
        }
    }
}
