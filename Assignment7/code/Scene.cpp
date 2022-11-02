//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;
    emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break;
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here

    Intersection inter_p = intersect(ray);
    // check if eye ray intersect with obj
    if(!inter_p.happened) 
        return {};
    // check if the intersection is a light source 
    if(inter_p.m->hasEmission()) 
        return inter_p.m->getEmission();

    Intersection inter;
    float pdf_light;
    sampleLight(inter, pdf_light);

    Vector3f L_dir, L_indir = 0.0;
    
    Vector3f wo(inter.coords - inter_p.coords); //light side

    Ray out_ray(inter_p.coords, wo.normalized());
    if(intersect(out_ray).distance > wo.norm() - 0.005f){
        L_dir = inter.emit * inter_p.m->eval(ray.direction, wo.normalized(), inter_p.normal) 
            * dotProduct(wo.normalized(), inter_p.normal)
            * dotProduct(-wo.normalized(), inter.normal)
            / std::pow(wo.norm(), 2) / pdf_light;
    }

    if(get_random_float() > RussianRoulette){
        return L_dir;
    }

    wo = inter_p.m->sample(ray.direction, inter_p.normal);
    out_ray = Ray(inter_p.coords, wo.normalized());

    // If out_ray hit a non-emitting object at q
    Intersection q = intersect(out_ray);
    float pdf = inter_p.m->pdf(ray.direction, wo.normalized(), inter_p.normal);
    if(q.happened && !q.m->hasEmission() && pdf > EPSILON){
        L_indir = castRay(out_ray, depth+1) * inter_p.m->eval(ray.direction, wo.normalized(), inter_p.normal)
            * dotProduct(wo.normalized(), inter_p.normal)
            / pdf
            / RussianRoulette;
    }

    return L_dir + L_indir;
}