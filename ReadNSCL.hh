#ifndef READNSCL_H
#define READNSCL_H
#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <string>
#include "stdint.h"

// Definitions
#define ALIGN8(x)    (((x)+7) & ~7)

// state change item type codes:
static const uint32_t BEGIN_RUN  = 1;
static const uint32_t END_RUN    = 2;
static const uint32_t PAUSE_RUN  = 3;
static const uint32_t RESUME_RUN = 4;

// Documentation item type codes:
static const uint32_t PACKET_TYPES        = 10;
static const uint32_t MONITORED_VARIABLES = 11;

// Scaler data:
static const uint32_t INCREMENTAL_SCALERS = 20;
static const uint32_t TIMESTAMPED_NONINCR_SCALERS =21;

// Physics events:
static const uint32_t PHYSICS_EVENT       = 30;
static const uint32_t PHYSICS_EVENT_COUNT = 31;

// Event builder related items:
static const uint32_t EVB_FRAGMENT        = 40; /* Event builder fragment. */
static const uint32_t EVB_UNKNOWN_PAYLOAD = 41; /* Evb fragment whose payload isn't a ring item */

// Classes
class EvtHeader
{
public:
  short int ID;   // 2 bytes
  short int Mask;
  int Serial;     // 4 bytes
  int Time;
  int Size;
};

class Bank
{
public:
  class BankHeader
  {
  public:
    char Name[4];
    int Type;
    int Size;
  };
  BankHeader header;
  std::vector<uint32_t> Data;

  void PrintBankHeader();
  void PrintBank();

};


// The main "event" class. Each event will have:
// A header and data
class Event {

public:

  int ReadEvent(std::ifstream& file);

  // Getters
  uint32_t getEventType(){return TypeCode;};
  uint32_t getEventSize(){return Size;};
  //  std::vector<uint16_t> getData(){return Data;};
  std::vector<Bank> getBanks(){return Banks;};

  // Setters

  // Printers
  void PrintEventHeader();
  void PrintEventData();


private:

  EvtHeader evtHead;
  std::vector<Bank> Banks;

  // Stored variables
  uint32_t Size;
  uint32_t TypeCode;
  std::vector<uint16_t> rawData;
  std::vector<uint16_t> Data;
 
  // Printers
  std::string TranslateType();

  
};

void ReadEventFile(const char* filename, std::vector<Event> *allEvents, int nEvents);

// define a collection of back and front position vectors
// struct _both_pos
// typedef struct _both_pos both_pos;


#endif
