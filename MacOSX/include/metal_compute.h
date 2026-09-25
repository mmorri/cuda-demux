#pragma once

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

namespace metal_compute {

class MetalBuffer {
public:
    id<MTLBuffer> buffer;
    size_t size;

    MetalBuffer(id<MTLDevice> device, size_t size_bytes, const void* data = nullptr) {
        size = size_bytes;
        if (data) {
            buffer = [device newBufferWithBytes:data length:size_bytes options:MTLResourceStorageModeShared];
        } else {
            buffer = [device newBufferWithLength:size_bytes options:MTLResourceStorageModeShared];
        }
    }

    void* contents() { return [buffer contents]; }
    const void* contents() const { return [buffer contents]; }
};

class MetalCompute {
private:
    id<MTLDevice> device;
    id<MTLCommandQueue> commandQueue;
    id<MTLLibrary> library;
    std::unordered_map<std::string, id<MTLComputePipelineState>> pipelines;

public:
    MetalCompute() {
        device = MTLCreateSystemDefaultDevice();
        if (!device) {
            throw std::runtime_error("Metal is not supported on this device");
        }
        commandQueue = [device newCommandQueue];
    }

    void loadShaderLibrary(const std::string& metallib_path) {
        NSString* path = [NSString stringWithUTF8String:metallib_path.c_str()];
        NSError* error = nil;
        library = [device newLibraryWithFile:path error:&error];

        if (!library) {
            NSString* errorStr = [error localizedDescription];
            throw std::runtime_error("Failed to load shader library: " + std::string([errorStr UTF8String]));
        }
    }

    void loadShaderSource(const std::string& source) {
        NSString* src = [NSString stringWithUTF8String:source.c_str()];
        NSError* error = nil;
        library = [device newLibraryWithSource:src options:nil error:&error];

        if (!library) {
            NSString* errorStr = [error localizedDescription];
            throw std::runtime_error("Failed to compile shader: " + std::string([errorStr UTF8String]));
        }
    }

    void createComputePipeline(const std::string& function_name) {
        NSString* funcName = [NSString stringWithUTF8String:function_name.c_str()];
        id<MTLFunction> function = [library newFunctionWithName:funcName];

        if (!function) {
            throw std::runtime_error("Failed to find function: " + function_name);
        }

        NSError* error = nil;
        id<MTLComputePipelineState> pipeline = [device newComputePipelineStateWithFunction:function error:&error];

        if (!pipeline) {
            NSString* errorStr = [error localizedDescription];
            throw std::runtime_error("Failed to create pipeline: " + std::string([errorStr UTF8String]));
        }

        pipelines[function_name] = pipeline;
    }

    std::shared_ptr<MetalBuffer> createBuffer(size_t size_bytes, const void* data = nullptr) {
        return std::make_shared<MetalBuffer>(device, size_bytes, data);
    }

    void dispatch(const std::string& kernel_name,
                  const std::vector<std::shared_ptr<MetalBuffer>>& buffers,
                  const std::vector<uint32_t>& constants,
                  MTLSize gridSize,
                  MTLSize threadGroupSize) {

        if (pipelines.find(kernel_name) == pipelines.end()) {
            throw std::runtime_error("Pipeline not found: " + kernel_name);
        }

        id<MTLComputePipelineState> pipeline = pipelines[kernel_name];
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];

        [encoder setComputePipelineState:pipeline];

        // Set buffers
        for (size_t i = 0; i < buffers.size(); i++) {
            [encoder setBuffer:buffers[i]->buffer offset:0 atIndex:i];
        }

        // Set constants
        size_t constant_offset = buffers.size();
        for (size_t i = 0; i < constants.size(); i++) {
            [encoder setBytes:&constants[i] length:sizeof(uint32_t) atIndex:constant_offset + i];
        }

        [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadGroupSize];
        [encoder endEncoding];
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];
    }

    void dispatchAsync(const std::string& kernel_name,
                       const std::vector<std::shared_ptr<MetalBuffer>>& buffers,
                       const std::vector<uint32_t>& constants,
                       MTLSize gridSize,
                       MTLSize threadGroupSize,
                       std::function<void()> completion) {

        if (pipelines.find(kernel_name) == pipelines.end()) {
            throw std::runtime_error("Pipeline not found: " + kernel_name);
        }

        id<MTLComputePipelineState> pipeline = pipelines[kernel_name];
        id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];

        [encoder setComputePipelineState:pipeline];

        // Set buffers
        for (size_t i = 0; i < buffers.size(); i++) {
            [encoder setBuffer:buffers[i]->buffer offset:0 atIndex:i];
        }

        // Set constants
        size_t constant_offset = buffers.size();
        for (size_t i = 0; i < constants.size(); i++) {
            [encoder setBytes:&constants[i] length:sizeof(uint32_t) atIndex:constant_offset + i];
        }

        [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadGroupSize];
        [encoder endEncoding];

        [commandBuffer addCompletedHandler:^(id<MTLCommandBuffer> buffer) {
            completion();
        }];

        [commandBuffer commit];
    }

    size_t getMaxThreadgroupSize(const std::string& kernel_name) {
        if (pipelines.find(kernel_name) == pipelines.end()) {
            throw std::runtime_error("Pipeline not found: " + kernel_name);
        }
        return [pipelines[kernel_name] maxTotalThreadsPerThreadgroup];
    }

    size_t getDeviceMemory() {
        return [device recommendedMaxWorkingSetSize];
    }
};

} // namespace metal_compute